/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <libbio/bed_reader.hh>
#include <libbio/file_handling.hh>
#include <libbio/utility/misc.hh> // log_time
#include <panvc3/msa_index.hh>
#include "cmdline.h"

namespace lb	= libbio;
namespace fs	= std::filesystem;


namespace {
	template <typename t_string>
	panvc3::msa_index::chr_entry const find_chr_entry(panvc3::msa_index::chr_entry_vector const &entries, t_string const &chr_id)
	{
		panvc3::msa_index::chr_entry_cmp chr_cmp;
		auto const it(std::lower_bound(entries.begin(), entries.end(), chr_id, chr_cmp));
		if (entries.end() == it || it->chr_id != chr_id)
		{
			std::cerr << "ERROR: Did not find an entry for chromosome ID “" << chr_id << "” in the MSA index." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		return *it;
	}
				
	
	template <typename t_string>
	panvc3::msa_index::sequence_entry const &find_sequence_entry(panvc3::msa_index::sequence_entry_vector const &entries, t_string const &seq_id)
	{
		panvc3::msa_index::sequence_entry_cmp seq_cmp;
		auto const it(std::lower_bound(entries.begin(), entries.end(), seq_id, seq_cmp));
		if (entries.end() == it || it->seq_id != seq_id)
		{
			std::cerr << "ERROR: Did not find an entry for sequence ID “" << seq_id << "” in the MSA index." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		return *it;
	}
	
	
	class bed_processor final : public lb::bed_reader_delegate
	{
	protected:
		panvc3::msa_index						m_msa_index;
		panvc3::msa_index::sequence_entry const	*m_src_entry{};
		panvc3::msa_index::sequence_entry const	*m_dst_entry{};
		std::string								m_chr_id;
		std::string								m_dst_chr_id;
		std::size_t								m_chr_id_matches{};
		std::size_t								m_chr_id_mismatches{};
		
	public:
		std::size_t convert_position(std::size_t const srcpos)
		{
			// Find the srcpos-th zero in the source sequence.
			auto const alnpos(m_src_entry->gap_positions_select0_support(1 + srcpos));
			
			// Check the dst character.
			bool const dstc(m_dst_entry->gap_positions[alnpos]);
			
			// Count the zeros in dst up to and including alnpos.
			// If dstc is zero, take the next character. Finally convert to zero-based.
			auto const dstpos(m_dst_entry->gap_positions_rank0_support(1 + alnpos) + dstc - 1);
			
			return dstpos;
		}
		
		
		void bed_reader_reported_error(std::size_t const lineno) override
		{
			std::cerr << "ERROR: Parse error in BED input on line " << lineno << '.' << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		
		// Reports a half-open interval.
		void bed_reader_found_region(std::string_view const chr_id, std::size_t const lb, std::size_t const rb) override
		{
			if (chr_id != m_chr_id)
			{
				++m_chr_id_mismatches;
				return;
			}
			
			++m_chr_id_matches;
			
			// Convert and output.
			auto const lb_(convert_position(lb));
			auto const rb_(convert_position(rb));
			std::cout << chr_id << '\t' << lb_ << '\t' << rb_ << '\n';
		}
		
		
		void bed_reader_did_finish() override
		{
			std::cout << std::flush;
			lb::log_time(std::cerr) << "Done. Chromosome ID matches: " << m_chr_id_matches << " mismatches: " << m_chr_id_mismatches << ".\n";
		}
		
		
		void process(std::istream &bed_stream, gengetopt_args_info &args_info)
		{
			m_chr_id = args_info.chr_arg;
			m_dst_chr_id = args_info.dst_chr_arg ?: "";
			
			lb::log_time(std::cerr) << "Loading the MSA index…\n";
		
			{
				lb::file_istream stream;
				lb::open_file_for_reading(args_info.msa_index_arg, stream);
				cereal::PortableBinaryInputArchive archive(stream);
				archive(m_msa_index);
			}
			
			// Find the relevant entries in the MSA index.
			{
				auto const &chr_entries(m_msa_index.chr_entries);
				auto const &src_chr(find_chr_entry(chr_entries, m_chr_id));
				auto const &dst_chr(m_dst_chr_id.empty() ? src_chr : find_chr_entry(chr_entries, m_dst_chr_id));
				auto const &src_seq_entries(src_chr.sequence_entries);
				auto const &dst_seq_entries(dst_chr.sequence_entries);
				m_src_entry = &find_sequence_entry(src_seq_entries, args_info.src_seq_arg);
				m_dst_entry = &find_sequence_entry(dst_seq_entries, args_info.dst_seq_arg);
			}
			
			// Process.
			lb::log_time(std::cerr) << "Processing the BED input…\n";
			lb::bed_reader bed_reader;
			bed_reader.read_regions(bed_stream, *this);
		}
	};
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	
	bed_processor processor;

	if (args_info.bed_arg)
	{
		lb::file_istream stream;
		lb::open_file_for_reading(args_info.bed_arg, stream);
		processor.process(stream, args_info);
	}
	else
	{
		processor.process(std::cin, args_info);
	}
	
	// Not reached.
	return EXIT_SUCCESS;
}
