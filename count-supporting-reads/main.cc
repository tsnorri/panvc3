/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <iostream>
#include <libbio/algorithm.hh>
#include <libbio/bed_reader.hh>
#include <libbio/vcf/region_variant_validator.hh>
#include <libbio/vcf/variant_end_pos.hh>
#include <libbio/vcf/vcf_reader.hh>
#include <libbio/file_handling.hh>
#include <panvc3/cigar.hh>
#include <panvc3/dna11_alphabet.hh>
#include <range/v3/all.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include "cmdline.h"

namespace fs	= std::filesystem;
namespace lb	= libbio;
namespace rsv	= ranges::views;
namespace vcf	= libbio::vcf;


namespace {

	typedef vcf::info_field <vcf::metadata_value_type::FLAG>	info_field_co;
	typedef vcf::info_field <vcf::metadata_value_type::FLAG>	info_field_usra;
	
	void output_cigar(std::vector <seqan3::cigar> const &cigar_seq, std::ostream &os)
	{
		using seqan3::get;
		
		for (auto const cigar_item : cigar_seq)
		{
			auto const op_count(get <0>(cigar_item));
			auto const operation(get <1>(cigar_item));
			os << op_count << operation.to_char();
		}
	}
	
	
	// From libvcf2multialign.
	// FIXME: come up with a way not to duplicate the code needed for storing field pointers.
	struct variant_format final : public vcf::variant_format
	{
		vcf::genotype_field_gt	*gt_field{};
		
		// Return a new empty instance of this class.
		virtual variant_format *new_instance() const override { return new variant_format(); }
		
		virtual void reader_did_update_format(vcf::reader &reader) override
		{
			this->assign_field_ptr("GT", gt_field);
		}
	};
	
	inline variant_format const &get_variant_format(vcf::variant const &var)
	{
		libbio_assert(var.reader()->has_assigned_variant_format());
		return static_cast <variant_format const &>(var.get_format());
	}
	
	inline variant_format const &get_variant_format(vcf::transient_variant const &var)
	{
		libbio_assert(var.reader()->has_assigned_variant_format());
		return static_cast <variant_format const &>(var.get_format());
	}
	
	
	class bed_reader_delegate : public vcf::region_variant_validator_bed_reader_delegate
	{
		typedef vcf::region_variant_validator_bed_reader_delegate base_t;
		
	public:
		using base_t::base_t;
		
		void bed_reader_reported_error(std::size_t const lineno) override
		{
			std::cerr << "ERROR: Parse error in BED on line " << lineno << ".\n";
			std::exit(EXIT_FAILURE);
		}
	};
	
	
	struct record_length
	{
		std::size_t reference_length{};
		std::size_t right_anchored_length{};
		
		record_length() = default;
		
		record_length(std::size_t const reference_length_, std::size_t const right_anchored_length_):
			reference_length(reference_length_),
			right_anchored_length(right_anchored_length_)
		{
		}
	};
	
	
	template <typename t_aln_record>
	record_length calculate_record_lengths(t_aln_record const &aln_record)
	{
		std::size_t reference_length{};
		std::size_t right_anchored_length{};
		//std::size_t query_length{};
		
		for (auto const &cigar_item : aln_record.cigar_sequence())
		{
			using seqan3::get;
			
			auto const op_count(get <0>(cigar_item));
			auto const operation(get <1>(cigar_item));
			auto const cc(operation.to_char());
			
			switch (cc)
			{
				case 'M': // Consume both query and reference.
				case '=':
				case 'X':
					reference_length += op_count;
					//query_length += op_count;
					right_anchored_length = reference_length;
					break;
				
				case 'I': // Consume query.
				case 'S':
					//query_length += op_count;
					break;
					
				case 'D': // Consume reference.
				case 'N':
					reference_length += op_count;
					break;
				
				case 'H':
				case 'P':
					break;
				
				default:
					libbio_fail("Unexpected CIGAR operation “", cc, "”");
					break;
			}
		}
		
		return {reference_length, right_anchored_length};
	}
	
	
	template <typename t_src, typename t_dst>
	inline void copy_sequence_slice(t_src const &src, std::size_t const pos, std::size_t const length, t_dst &dst)
	{
		ranges::copy(src | rsv::slice(pos, pos + length), ranges::back_inserter(dst));
	}
	
	
	template <typename t_src, typename t_dst>
	inline void copy_sequence_slice_clipped(t_src const &src, std::size_t const pos, std::size_t const length, t_dst &dst)
	{
		// Construct clipped panvc3::dna11 characters from plain seqan3::dna5 characters.
		ranges::copy(
			src
			| rsv::slice(pos, pos + length)
			| rsv::transform([](auto const &cc){ return panvc3::dna11(cc, true); }),
			ranges::back_inserter(dst)
		);
	}


	// Check compatibility in terms of what operations can be used to continue
	bool can_continue_with_operation(seqan3::cigar::operation const lhs, seqan3::cigar::operation const rhs)
	{
		// The specification forbids some combinations, e.g. H can only appear as the first
		// or the last operation. We handle them anyway b.c. it is not very complicated.
		switch (lhs.to_char())
		{
			case 'D':
			case 'N':
			{
				auto const cc(rhs.to_char());
				return ('D' == cc || 'N' == cc);
			}
				
			case 'H':
			case 'P':
			{
				auto const cc(rhs.to_char());
				return ('H' == cc || 'P' == cc);
			}
				
			default:
				break;
		}
		
		return lhs == rhs;
	}

	
	template <typename t_aln_record> // FIXME: try to find a way to get the sequence type crom t_aln_record instead of assuming dna5.
	bool try_read_aligned_sequence(
		t_aln_record const &aln_record,
		std::size_t const rec_len,
		std::size_t const var_pos,
		std::size_t var_ref_len,
		std::size_t var_alt_len,
		panvc3::dna11_vector &dst,
		bool const should_include_clipping
	)
	{
		// Caller checks that the variant is contained in the read.
		using panvc3::operator""_dna11;
		using seqan3::literals::operator""_cigar_operation;
		
		dst.clear();
		auto const var_ref_len_(var_ref_len);
		auto const var_alt_len_(var_alt_len);
		
		// Process the CIGAR string and read the characters the alignment suggests.
		auto const &seq(aln_record.sequence()); // std::vector <seqan3::dna5> by default.
		auto const &cigar_seq(aln_record.cigar_sequence());
		auto it(cigar_seq.begin());
		auto const end(cigar_seq.end());

		try
		{
			auto const rec_pos_(aln_record.reference_position());
			libbio_assert(rec_pos_);
			libbio_assert_lte(0, *rec_pos_);
			std::size_t rec_pos(*rec_pos_);
			libbio_assert_lte(rec_pos, var_pos);
			libbio_assert_lte(var_pos + var_ref_len, rec_pos + rec_len);
			
			panvc3::cigar_count_type op_count{};
			seqan3::cigar::operation operation{};
			std::size_t seg_pos{};
			for (; it != end; ++it)
			{
				using seqan3::get;
				op_count = get <0>(*it);
				operation = get <1>(*it);
				
				if (rec_pos < var_pos)
				{
					// Find the variant position in the aligned segment.
					auto const operation_(operation.to_char());
					switch (operation_)
					{
						case 'M':	// Consume both query and reference.
						case '=':
						case 'X':
						{
							auto const min_length(lb::min_ct(var_pos - rec_pos, op_count));
							op_count -= min_length;
							rec_pos += min_length;
							seg_pos += min_length;
							if (op_count)
								break;
							
							continue; // Continues the enclosing loop.
						}
					
						case 'D':	// Consume reference.
						case 'N':
						{
							auto const min_length(lb::min_ct(var_pos - rec_pos, op_count));
							op_count -= min_length;
							rec_pos += min_length;
							if (op_count)
								break;
						
							continue; // Continues the enclosing loop.
						}
						
						case 'I':	// Consume query.
						case 'S':
						{
							// Not yet at var_pos, so we need to skip this segment.
							seg_pos += op_count;
							continue; // Continues the enclosing loop.
						}
						
						default:
							continue; // Continues the enclosing loop.
					}
				}
				
				// Read the aligned sequence.
				// Reached iff. var_pos ≤ rec_pos or op_count is non-zero as a result of the last subtraction.
				auto const operation_(operation.to_char());
				switch (operation_)
				{
					case 'M': // Consume both query and reference.
					case '=':
					case 'X':
					{
						// Since the segment is aligned to the reference, consume as many
						// characters that correspond to REF as possible.
						auto const min_length(lb::min_ct(op_count, var_ref_len));
						copy_sequence_slice(seq, seg_pos, min_length, dst);
						seg_pos += min_length;
						var_ref_len -= std::min(var_ref_len, min_length);
						var_alt_len -= std::min(var_alt_len, min_length);
						op_count -= min_length;
						
						// Handle the case where the variant has an insertion (relative to the
						// current position) but the whole M / = / X operation was not consumed.
						if (op_count && 0 == var_ref_len)
							var_alt_len = 0; // Clearly the read does not have an insertion since the nucleotides are aligned.
						
						break;
					}
					
					case 'I': // Consumes query.
					{
						// Consider the whole insertion.
						copy_sequence_slice(seq, seg_pos, op_count, dst);
						seg_pos += op_count;
						var_alt_len -= lb::min_ct(var_alt_len, op_count);
						op_count = 0;
						break;
					}
					
					case 'D': // Consume reference.
					case 'N':
					{
						if (var_ref_len < op_count && (dst.empty() || dst.back() != '~'_dna11))
							dst.push_back('~'_dna11);
						
						var_ref_len -= lb::min_ct(var_ref_len, op_count);
						op_count = 0;
						break;
					}
					
					case 'S': // Consume query.
					{
						if (should_include_clipping)
							copy_sequence_slice_clipped(seq, seg_pos, op_count, dst);
						seg_pos += op_count;
						var_alt_len -= lb::min_ct(var_alt_len, op_count); // FIXME: Should the characters actually be considered part of ALT?
						op_count = 0;
						break;
					}
					
					case 'H': // Consume nothing.
					case 'P':
						op_count = 0;
						break;
						
					default:
						libbio_fail("Unexpected CIGAR operation “", operation_, "”");
						break;
				}
				
				if (0 == var_ref_len && 0 == var_alt_len)
					goto finish; // Skip incrementing the iterator.
			}
			
			// This should be the case now.
			libbio_assert(0 != var_ref_len || 0 != var_alt_len);
			return false;
			
		finish:
			libbio_assert(0 == var_ref_len && 0 == var_alt_len);

			// The SAM/BAM specification only states that adjacent CIGAR operations “should” be different
			// (Sequence Alignment/Map Format Specification, 22nd Aug 2022, Section 2, Item 2),
			// so we try to handle the trailing operations that match if needed.
			if (0 == op_count && it != end) // op_count ≠ 0 => There are remaining M / = / X operations that need not be handled.
			{
				// Check the next operation if there is one.
				++it;
				auto const prev_operation(operation);
				for (; it != end; ++it)
				{
					using seqan3::get;
					op_count = get <0>(*it);
					operation = get <1>(*it);
					
					if ('S'_cigar_operation == operation)
					{
						if (should_include_clipping)
							copy_sequence_slice_clipped(seq, seg_pos, op_count, dst);
						seg_pos += op_count;
						continue;
					}
					
					if (!can_continue_with_operation(prev_operation, operation))
						break;
					
					// Continue with the same operation.
					auto const operation_(operation.to_char());
					switch (operation_)
					{
						case 'I': // Consumes query.
						{
							copy_sequence_slice(seq, seg_pos, op_count, dst);
							seg_pos += op_count;
							break;
						}
						
						case 'D': // Consume reference.
						case 'N':
						{
							if (dst.empty() || dst.back() != '~'_dna11)
								dst.push_back('~'_dna11);
							break;
						}
						
						case 'H': // Consume nothing.
						case 'P':
							break;
							
						default:
							libbio_fail("Unexpected CIGAR operation “", operation_, "”");
							break;
					}
				}
			}
			
			return true;
		}
		catch (lb::assertion_failure_exception const &exc)
		{
			// For debugging.
			std::cerr << "CIGAR: "; 
			output_cigar(cigar_seq, std::cerr);
			std::cerr << '\n';
			throw;
		}
	}
	
	
	// We use this structure to preprocess the given alignment record.
	// Currently the reference length is cached, since determining it from the alignment record
	// itself seems difficult.
	template <typename t_aln_record>
	struct record
	{
		t_aln_record	alignment_record;
		record_length	rec_length;
		
		record(t_aln_record &alignment_record_, record_length const rec_length_):
			alignment_record(alignment_record_),
			rec_length(rec_length_)
		{
			libbio_assert(alignment_record.reference_position().has_value()); // Make sure that the optional holds a value.
		}

		explicit record(t_aln_record &alignment_record_):
			record(alignment_record_, calculate_record_lengths(alignment_record))
		{
		}
		
		auto reference_position() const { return *alignment_record.reference_position(); }
		std::size_t reference_length() const { return rec_length.reference_length; }
		std::size_t reference_end() const { return reference_position() + reference_length(); }
		std::size_t right_anchored_length() const { return rec_length.right_anchored_length; }
		
		// FIXME: try to find a way to get the sequence type from t_aln_record, i.e. dna11 only if the source is dna5.
		bool try_read_aligned_sequence(
			std::size_t const var_pos,
			std::size_t const var_ref_len,
			std::size_t const var_alt_len,
			panvc3::dna11_vector &dst,
			bool const should_include_clipping
		) const
		{
			return ::try_read_aligned_sequence(alignment_record, reference_length(), var_pos, var_ref_len, var_alt_len, dst, should_include_clipping);
		}
	};
	
	
	struct reference_position_cmp
	{
		template <typename t_rec>
		bool operator()(t_rec const &lhs, t_rec const &rhs) const
		{
			return lhs.reference_position() < rhs.reference_position();
		}
	};
	
	
	constexpr static auto sam_reader_fields()
	{
		return seqan3::fields <
			seqan3::field::ref_id,
			seqan3::field::ref_offset,
			seqan3::field::mate,
			seqan3::field::seq,
			seqan3::field::flag,
			seqan3::field::cigar
		>{};
	}
	
	
	static auto open_alignment_input_file(fs::path const &path)
	{
		return seqan3::sam_file_input(path, sam_reader_fields());
	}
	
	
	template <typename t_format>
	static auto open_alignment_input_stream(std::istream &stream, t_format &&format)
	{
		return seqan3::sam_file_input(stream, std::forward <t_format>(format), sam_reader_fields());
	}
	
	
	template <typename t_aln_input>
	class alignment_reader
	{
	public:
		typedef typename t_aln_input::traits_type				alignment_input_traits_type;
		typedef typename alignment_input_traits_type::ref_ids	reference_ids_type;
		typedef typename t_aln_input::iterator					alignment_iterator_type;
		typedef typename t_aln_input::sentinel					alignment_sentinel_type;
		typedef typename t_aln_input::record_type				alignment_record_type;
		typedef record <alignment_record_type>					record_type;
		typedef std::set <record_type, reference_position_cmp>	record_set;
		
		struct alignment_statistics
		{
			std::size_t	reads_processed{};
			std::size_t	flags_not_matched{};
			std::size_t	ref_id_mismatches{};
			std::size_t	mate_ref_id_mismatches{};
			std::size_t	position_mismatches{};
			std::size_t	matched_reads{};
		};
		
	protected:
		alignment_iterator_type		m_it;
		alignment_sentinel_type		m_end;
		reference_ids_type const	&m_reference_ids; // Change to a pointer if needed.
		char const					*m_contig_prefix{};
		record_set					m_candidate_records;
		alignment_statistics		m_statistics{};
		std::size_t					m_prev_record_pos{};
		bool						m_should_consider_primary_alignments_only{};
		bool						m_requires_same_config_prefix_in_next{};
		
		
	public:
		alignment_reader(
			t_aln_input &aln_input,
			char const *contig_prefix,
			bool const should_consider_primary_alignments_only,
			bool const requires_same_config_prefix_in_next
		):
			m_it(aln_input.begin()),
			m_end(aln_input.end()),
			m_reference_ids(aln_input.header().ref_ids()),
			m_contig_prefix(contig_prefix),
			m_should_consider_primary_alignments_only(should_consider_primary_alignments_only),
			m_requires_same_config_prefix_in_next(requires_same_config_prefix_in_next)
		{
		}
		
		record_set const &candidate_records() const { return m_candidate_records; }
		alignment_statistics statistics() const { return m_statistics; }
		
		void update_candidate_records(std::size_t const var_pos)
		{
			// Remove the old records.
			std::erase_if(m_candidate_records, [var_pos](auto const &rec){
				return rec.reference_end() <= var_pos;
			});
			
			// Iterate over the alignment records.
			for (; m_it != m_end; ++m_it)
			{
				++m_statistics.reads_processed;
				if (0 == m_statistics.reads_processed % 10000000)
					lb::log_time(std::cerr) << "Processed " << m_statistics.reads_processed << " reads…\n";
				
				// If POS is zero, skip the record.
				auto &aln_rec(*m_it);
				auto const flags(aln_rec.flag());
				
				if (lb::to_underlying(flags & (
					seqan3::sam_flag::unmapped					|
					seqan3::sam_flag::failed_filter				|
					seqan3::sam_flag::duplicate					|
					seqan3::sam_flag::supplementary_alignment
				))) // Ignore unmapped, filtered, duplicate and supplementary.
				{
					++m_statistics.flags_not_matched;
					continue;
				}
				
				// Ignore secondary if requested.
				if (m_should_consider_primary_alignments_only && lb::to_underlying(flags & seqan3::sam_flag::secondary_alignment))
				{
					++m_statistics.flags_not_matched;
					continue;
				}
				
				// Check the contig name if requested.
				auto const ref_id(aln_rec.reference_id());
				if (m_contig_prefix)
				{
					// Check the contig prefix for the current alignment.
					if (!ref_id.has_value())
					{
						++m_statistics.ref_id_mismatches;
						continue;
					}
					
					if (!m_reference_ids[*ref_id].starts_with(m_contig_prefix))
					{
						++m_statistics.ref_id_mismatches;
						continue;
					}
					
					// Check the contig prefix of the next primary alignment.
					// (Currently we do not try to determine the reference ID of all the possible alignments
					// of the next read.)
					if (m_requires_same_config_prefix_in_next)
					{
						auto const mate_ref_id(aln_rec.mate_reference_id());
						if (!mate_ref_id.has_value())
						{
							++m_statistics.mate_ref_id_mismatches;
							continue;
						}
						
						if (!m_reference_ids[*mate_ref_id].starts_with(m_contig_prefix))
						{
							++m_statistics.mate_ref_id_mismatches;
							continue;
						}
					}
				}
				
				auto const rec_ref_pos_(aln_rec.reference_position());
				if (!rec_ref_pos_.has_value())
				{
					++m_statistics.flags_not_matched;
					continue;
				}
				
				auto const rec_ref_pos(*rec_ref_pos_);
				// The following assertion checks that rec_ref_pos is non-negative, too.
				libbio_always_assert_lte(m_prev_record_pos, rec_ref_pos); // FIXME: error message: alignments need to be sorted by position.
				
				m_prev_record_pos = rec_ref_pos; // We could also consider some of the filtered reads’ positions for m_prev_record_pos, but I think this will suffice.
				
				// Check that the alignment starts before the left boundary of the variant.
				if (! (std::size_t(rec_ref_pos) <= var_pos))
				{
					++m_statistics.position_mismatches;
					return;
				}
				
				auto const rec_lengths(calculate_record_lengths(aln_rec));
				auto const rec_ref_end(rec_ref_pos + rec_lengths.reference_length);
				if (rec_ref_end <= var_pos)
				{
					++m_statistics.position_mismatches;
					continue;
				}
				
				// We could calculate the number of matched reads from the other statistics
				// but having a separate variable is safer w.r.t. potential bugs.
				++m_statistics.matched_reads;
				m_candidate_records.emplace(aln_rec, rec_lengths);
			}
		}
	};
	
	
	struct variant_statistics
	{
		std::size_t						variants_processed{};
		std::size_t						chr_id_mismatches{};
		std::size_t						zygosity_mismatches{};
		std::size_t						zero_coverage{};
	};
	
	
	class variant_validator final : public vcf::region_variant_validator
	{
	public:
		vcf::variant_validation_result handle_unordered_contigs(vcf::transient_variant const &var) override
		{
			std::cerr << "ERROR: Line " << var.lineno() << ": Variants are not sorted by chromosome ID and position.\n";
			std::exit(EXIT_FAILURE);
		}
		
		vcf::variant_validation_result handle_unordered_variants(vcf::transient_variant const &var) override
		{
			std::cerr << "ERROR: Line " << var.lineno() << ": Contigs are not in contiguous blocks.\n";
			std::exit(EXIT_FAILURE);
		}
	};


	void output_alts(vcf::transient_variant const &var, std::ostream &os)
	{
		bool is_first(true);
		auto const &alts(var.alts());
		for (auto const &alt : alts)
		{
			if (!is_first)
				os << ',';
			
			switch (alt.alt_sv_type)
			{
				case vcf::sv_type::NONE:
					os << alt.alt;
					break;
					
				case vcf::sv_type::UNKNOWN:
					continue;
				
				case vcf::sv_type::DEL:
				case vcf::sv_type::DEL_ME:
					os << "<DEL>";
					break;
					
				case vcf::sv_type::INS:
				case vcf::sv_type::DUP:
				case vcf::sv_type::INV:
				case vcf::sv_type::CNV:
				case vcf::sv_type::DUP_TANDEM:
				case vcf::sv_type::INS_ME:
				case vcf::sv_type::UNKNOWN_SV:
					throw std::runtime_error("ALT type not handled");
			}
			
			is_first = false;
		}
	}
	
	
	template <typename t_aln_input>
	void process_(
		t_aln_input &aln_input,
		gengetopt_args_info const &args_info
	)
	{
		// Overview of the algorithm
		// =========================
		// – Process the VCF and SAM/BAM records from the smallest co-ordinate to the largest.
		// – For each VCF record, determine the set of alignments where the REF sequence of the variant
		//   is fully contained in the alignment. Process each such alignment.
		// – Find the starting point of the variant in question within the alignment.
		// – Copy characters from the read as indicated by the CIGAR operations until both |REF| and
		//   |ALT| have been reached. The only exception is that the process is stopped early if the
		//   next aligned pair of characters (i.e. an M / = / X operation) is found.
		// – If the rightmost operation was an insertion or a deletion, continue handling similar
		//   operations.
		// – If both limits were reached, the operation is consiered successful.
		
		auto const chr_id(args_info.chr_arg);
		auto const vcf_path(args_info.vcf_arg);
		auto const regions_path(args_info.regions_arg);
		bool const should_include_clipping(args_info.include_clipping_flag);
		bool const should_consider_primary_alignments_only(args_info.primary_only_flag);
		bool const requires_same_config_prefix_in_next(args_info.same_ref_flag);
		bool const should_anchor_reads_left_only(args_info.anchor_left_flag);
		
		// Open the VCF file.
		vcf::stream_input <lb::file_istream> vcf_input;
		lb::open_file_for_reading(vcf_path, vcf_input.stream());
		vcf::reader vcf_reader(vcf_input);
		
		vcf::add_reserved_info_keys(vcf_reader.info_fields());
		vcf::add_reserved_genotype_keys(vcf_reader.genotype_fields());
		
		// Read the headers.
		vcf_reader.set_variant_format(new variant_format());
		vcf_reader.read_header();
		
		// Retrieve the END field.
		vcf::info_field_end *vcf_end_field{};
		info_field_co *vcf_co_field{};
		info_field_usra *vcf_usra_field{};
		vcf_reader.get_info_field_ptr(args_info.end_field_id_arg, vcf_end_field);
		vcf_reader.get_info_field_ptr(args_info.co_field_id_arg, vcf_co_field);
		vcf_reader.get_info_field_ptr(args_info.usra_field_id_arg, vcf_usra_field);
		libbio_always_assert(vcf_end_field); // Reserved so should be set.
		
		// Read the regions if needed, otherwise make sure that the variants are in a consistent order.
		vcf::region_variant_validator variant_validator(nullptr != regions_path);
		vcf_reader.set_variant_validator(variant_validator);
		
		{
			lb::bed_reader bed_reader;
			bed_reader_delegate delegate(variant_validator.regions());
			bed_reader.read_regions(regions_path, delegate);
		}
		
		alignment_reader aln_reader(
			aln_input,
			args_info.contig_prefix_arg,
			should_consider_primary_alignments_only,
			requires_same_config_prefix_in_next
		);
		std::map <panvc3::dna11_vector, std::size_t> supported_sequences;
		panvc3::dna11_vector buffer;
		variant_statistics var_statistics;
		
		vcf_reader.set_parsed_fields(vcf::field::ALL);
		lb::log_time(std::cerr) << "Processing…\n";
		vcf_reader.parse(
			[
				chr_id,							// Pointer
				vcf_end_field,					// Pointer
				vcf_co_field,					// Pointer
				vcf_usra_field,					// Pointer
				should_include_clipping,		// bool
				should_anchor_reads_left_only,	// bool
				&aln_reader,
				&supported_sequences,
				&buffer,
				&var_statistics
			](vcf::transient_variant const &var){
				++var_statistics.variants_processed;
				if (0 == var_statistics.variants_processed % 100000)
					lb::log_time(std::cerr) << "Processed " << var_statistics.variants_processed << " variants…\n";
				
				auto const var_pos(var.zero_based_pos());
				
				auto const gt_field_ptr(get_variant_format(var).gt_field);
				libbio_always_assert(gt_field_ptr); // The variants should always have the GT field.
				auto const &gt_field(*gt_field_ptr);
				
				auto const &samples(var.samples());
				libbio_always_assert_eq(1, samples.size()); // FIXME: handle more than one sample, e.g by allowing the user to specify the sample identifier.
				
				// Check the chromosome identifier.
				if (chr_id && (chr_id != var.chrom_id()))
				{
					++var_statistics.chr_id_mismatches;
					return true;
				}
				
				// Get the sample genotype.
				auto const &sample(samples.front());
				auto const &gt(gt_field(sample)); // vector of sample_genotype
				libbio_always_assert_eq_msg(2, gt.size(), "Variant on line ", var.lineno(), " has non-diploid GT (", gt.size(), ")"); // FIXME: error message, or handle other zygosities.
				
				// Check the zygosity. (Generalized for polyploid.)
				static_assert(0x7fff == vcf::sample_genotype::NULL_ALLELE); // Should be positive and small enough s.t. the sum can fit into std::uint64_t or similar.
				auto const zygosity(std::accumulate(gt.begin(), gt.end(), std::uint64_t(0), [](auto const acc, vcf::sample_genotype const &sample_gt){
					return acc + (sample_gt.alt ? 1 : 0);
				}));
				
				if (1 != zygosity)
				{
					++var_statistics.zygosity_mismatches;
					return true;
				}
				
				// A suitable variant was found.
				auto const var_end_pos(vcf::variant_end_pos(var, *vcf_end_field));
				libbio_assert_lte(var_pos, var_end_pos);
				
				// Update the set of matching alignments.
				aln_reader.update_candidate_records(var_pos);
				auto const &candidate_records(aln_reader.candidate_records());
				if (candidate_records.empty())
				{
					++var_statistics.zero_coverage;
					return true;
				}
				
				// Output “V” chrom pos id(s) ref alts is_reversed, separated by tabs.
				std::cout << "V\t" << var.chrom_id() << '\t' << var_pos << '\t';
				ranges::copy(var.id(), ranges::make_ostream_joiner(std::cout, ","));
				std::cout << '\t' << var.ref() << '\t';
				output_alts(var, std::cout);
				std::cout << '\t' << +((vcf_co_field && vcf_co_field->has_value(var)) || (vcf_usra_field && vcf_usra_field->has_value(var))) << '\n';
				
				// FIXME: While we output all ALTs above, considering each for counting supporting reads would be quite difficult. Since we only handle heterozygous variants of diploid donors (for now), this is not needed.
				libbio_always_assert_eq(1, var.alts().size());
				auto const &alt(var.alts().front());
				auto const var_alt_len(alt.alt.size());
				
				auto const var_ref_len(var_end_pos - var_pos);
				
				supported_sequences.clear();
				for (auto const &rec : candidate_records)
				{
					auto const rec_pos(rec.reference_position());
					
					try
					{
						libbio_assert_lte(rec_pos, var_pos);
						
						// Add to the supported sequences if possible.
						// Check first that at least the reference sequence is present in the current read.
						// FIXME: Do we need to count the mismatch? By variant?
						if ((
								should_anchor_reads_left_only
								? var_end_pos <= rec_pos + rec.reference_length()
								: var_end_pos <  rec_pos + rec.right_anchored_length()
							) && rec.try_read_aligned_sequence(var_pos, var_ref_len, var_alt_len, buffer, should_include_clipping)
						)
						{
							auto &&[it, did_emplace] = supported_sequences.try_emplace(buffer, 0);
							++(it->second);
						}
					}
					catch (lb::assertion_failure_exception const &exc)
					{
						// For debugging.
						std::cerr << "rec_ref_len: " << rec.reference_length() << '\n';
						std::cerr << "rec_pos:     " << rec_pos << '\n';
						std::cerr << "var_pos:     " << var_pos << '\n';
						std::cerr << "var_end_pos: " << var_end_pos << '\n';
						std::cerr << "var_ref_len: " << var_ref_len << '\n';
						std::cerr << "var_alt_len: " << var_alt_len << '\n';
						std::cerr << "ALT:         " << alt.alt << '\n';
						throw;
					}
				}
				
				// Output “R” coverage sequence, separated by tabs for each distinct subsequence.
				for (auto const &kv : supported_sequences)
				{
					std::cout << "R\t" << kv.second << '\t';
					if (kv.first.empty())
						std::cout << "<DEL>";
					else
						std::copy(kv.first.begin(), kv.first.end(), std::ostream_iterator <panvc3::dna11>(std::cout));
					std::cout << '\n';
				}
				
				return true;
			}
		);
		
		// Output statistics.
		{
			auto const chr_id_mismatches(var_statistics.chr_id_mismatches + variant_validator.chromosome_id_mismatches());
			
			std::cout << "S\tTotal variants\t"				<< var_statistics.variants_processed		<< '\n';
			std::cout << "S\tChromosome ID mismatches\t"	<< chr_id_mismatches						<< '\n';
			std::cout << "S\tPosition mismatches\t"			<< variant_validator.position_mismatches()	<< '\n';
			std::cout << "S\tZygosity mismatches\t"			<< var_statistics.zygosity_mismatches		<< '\n';
			std::cout << "S\tZero coverage\t"				<< var_statistics.zero_coverage				<< '\n';
			std::cout << std::flush;
		}
		
		{
			auto const &aln_statistics(aln_reader.statistics());
			std::cout << "T\tReads processed\t"				<< aln_statistics.reads_processed			<< '\n';
			std::cout << "T\tFlags not matched\t"			<< aln_statistics.flags_not_matched			<< '\n';
			std::cout << "T\tRef. ID mismatches\t"			<< aln_statistics.ref_id_mismatches			<< '\n';
			std::cout << "T\tPair ref. ID mismatches\t"		<< aln_statistics.mate_ref_id_mismatches	<< '\n';
			std::cout << "T\tPosition mismatches\t"			<< aln_statistics.position_mismatches		<< '\n';
			std::cout << "T\tMatched alignments\t"			<< aln_statistics.matched_reads				<< '\n';
		}
		
		lb::log_time(std::cerr) << "Done.\n";
	}


	void process(gengetopt_args_info const &args_info)
	{
		// Open the SAM input. We expect the alignements to have been sorted by the leftmost co-ordinate.
		if (args_info.alignments_arg)
		{
			fs::path const alignments_path(args_info.alignments_arg);
			auto aln_input(open_alignment_input_file(alignments_path));
			process_(aln_input, args_info);
		}
		else
		{
			auto aln_input(open_alignment_input_stream(std::cin, seqan3::format_sam{}));
			process_(aln_input, args_info);
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	
	if (args_info.same_ref_flag)
	{
		if (!args_info.contig_prefix_arg)
		{
			std::cerr << "ERROR: --same-ref requires --contig-prefix." << std::endl;
			return EXIT_FAILURE;
		}
		
		if (!args_info.primary_only_flag)
		{
			std::cerr << "ERROR: --same-ref currently requires --primary-only." << std::endl;
			return EXIT_FAILURE;
		}
	}
	
	process(args_info);
	
	return EXIT_SUCCESS;
}
