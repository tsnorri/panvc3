/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/dispatch.hh>
#include <panvc3/msa_index.hh>
#include "index_handling.hh"
#include "input_processor.hh"

namespace lb	= libbio;


namespace {
	
	// Helpers for the sequence entry list input.
	struct sequence_entry
	{
		std::string	seq_id;
		std::string	path;

		sequence_entry(std::string const &seq_id_, std::string const &path_):
			seq_id(seq_id_),
			path(path_)
		{
		}
	};
	
	
	struct chr_entry
	{
		std::string						chr_id;
		std::vector <sequence_entry>	sequence_entries;
		
		explicit chr_entry(std::string const &chr_id_):
			chr_id(chr_id_)
		{
		}
		
		bool operator<(chr_entry const &other) const { return chr_id < other.chr_id; }
	};
	
	typedef std::vector <chr_entry> chr_entry_vector;
	
	struct chr_entry_cmp
	{
		bool operator()(chr_entry const &lhs, std::string const &rhs) const { return lhs.chr_id < rhs; }
		bool operator()(std::string const &lhs, chr_entry const &rhs) const { return lhs < rhs.chr_id; }
	};
	
	
	template <typename t_type, typename t_cmp>
	t_type &find_or_insert(std::vector <t_type> &vec, std::string const &key, t_cmp const &cmp)
	{
		auto const rng(std::equal_range(vec.begin(), vec.end(), key, cmp));
		if (rng.first == rng.second)
			return *vec.emplace(rng.second, key);
		else
			return *rng.first;
	}
	
	
	panvc3::msa_index::chr_entry &find_msa_chr_entry(panvc3::msa_index &msa_index, std::string const &chr_id)
	{
		panvc3::msa_index::chr_entry_cmp cmp;
		return find_or_insert(msa_index.chr_entries, chr_id, cmp);
	}
	
	
	chr_entry &find_chr_entry(chr_entry_vector &vec, std::string const &chr_id)
	{
		chr_entry_cmp cmp;
		return find_or_insert(vec, chr_id, cmp);
	}
	
	
	void read_input_entry(std::size_t const lineno, std::string const &entry, std::string &chr_id, std::string &seq_id, std::string &path)
	{
		std::string_view const &entry_sv(entry);
		std::size_t start(0);

		std::array dsts{&chr_id, &seq_id};
		
		for (auto *dst_ptr : dsts)
		{
			auto const tab_pos(entry_sv.find('\t', start));
			if (std::string::npos == tab_pos)
			{
				std::cerr << "ERROR: Parse error in input on line " << lineno << ".\n";
				std::exit(EXIT_FAILURE);
			}
			
			*dst_ptr = entry_sv.substr(start, tab_pos - start);
			start = 1 + tab_pos;
		}
		
		path = entry_sv.substr(start);
	}
}


namespace panvc3::msa_indices {
	
	void sequence_list_input_processor::process(input_handler &handler)
	{
		panvc3::msa_index msa_index;
			
		lb::file_istream path_stream;
		lb::file_ostream msa_index_stream;

		if (!m_msa_index_input_path.empty())
			load_msa_index(m_msa_index_input_path.c_str(), msa_index);
		
		lb::log_time(std::cerr) << "Loading the input sequences…\n";
		lb::open_file_for_reading(m_input_path, path_stream);
		lb::open_file_for_writing(m_msa_index_output_path, msa_index_stream, lb::make_writing_open_mode({lb::writing_open_mode::CREATE})); // FIXME: add overwriting conditionally.
		cereal::PortableBinaryOutputArchive msa_archive(msa_index_stream);
		
		chr_entry_vector chr_entries;
		
		// Read the paths.
		{
			std::string entry_buffer;
			std::string chr_id;
			std::string seq_id;
			std::string path;
			std::size_t lineno(1);
			
			while (std::getline(path_stream, entry_buffer))
			{
				read_input_entry(lineno, entry_buffer, chr_id, seq_id, path);
				auto &entry(find_chr_entry(chr_entries, chr_id));
				entry.sequence_entries.emplace_back(seq_id, path);
				++lineno;
			}
		}
		
		// Prepare the MSA index.
		msa_index.chr_entries.reserve(msa_index.chr_entries.size() + chr_entries.size());
		
		// Handle the inputs and compress and sort in background.
		std::vector <char> seq_buffer;
		lb::dispatch_ptr <dispatch_group_t> main_group(dispatch_group_create());
		auto *global_queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));

		for (auto const &chr_entry : chr_entries)
		{
			lb::log_time(std::cerr) << "Handling sequences for chromosome " << chr_entry.chr_id << "…\n";
			
			lb::dispatch_ptr <dispatch_group_t> chr_group(dispatch_group_create());
			auto &msa_chr_entry(find_msa_chr_entry(msa_index, chr_entry.chr_id));
			msa_chr_entry.sequence_entries.reserve(chr_entry.sequence_entries.size());
			
			std::size_t input_size{};
			for (auto const &seq_entry : chr_entry.sequence_entries)
			{
				handler.process_input(
					seq_entry.path,
					[
						this,
						global_queue,
						&seq_buffer,
						&chr_entry,
						&seq_entry,
						&chr_group,
						&msa_chr_entry,
						&input_size
					](lb::file_handle &handle){
						sdsl::bit_vector bv;
						lb::log_time(std::cerr) << "Processing " << seq_entry.path << "…\n";
						build_index_vector_one_sequence(m_index_vector_builder, chr_entry.chr_id, seq_entry.seq_id, seq_buffer, handle, input_size);
						
						auto &msa_seq_entry(msa_chr_entry.sequence_entries.emplace_back());
						lb::dispatch_group_async_fn(*chr_group, global_queue, [
							bv = m_index_vector_builder.destination_bit_vector(),
							&seq_entry,
							&msa_seq_entry
						](){
							msa_seq_entry = panvc3::msa_index::sequence_entry(seq_entry.seq_id, bv);
						});
					}
				);
			}
			
			dispatch_group_enter(*main_group);
			{
				auto main_group_(*main_group);
				dispatch_group_notify(*chr_group, global_queue, ^{
					std::sort(msa_chr_entry.sequence_entries.begin(), msa_chr_entry.sequence_entries.end());
					dispatch_group_leave(main_group_);
				});
			}
		}
		
		lb::log_time(std::cerr) << "Compressing the MSA index…\n";
		dispatch_group_wait(*main_group, DISPATCH_TIME_FOREVER);
		lb::log_time(std::cerr) << "Sorting the remaining index entries…\n";
		std::sort(msa_index.chr_entries.begin(), msa_index.chr_entries.end());
		lb::log_time(std::cerr) << "Serialising the MSA index…\n";
		msa_archive(msa_index);
		msa_index_stream << std::flush;
		lb::log_time(std::cerr) << "Done.\n";
		std::exit(EXIT_SUCCESS);
	}
	
	
	void a2m_input_processor::index_vector_builder_did_process_sequence(
		index_vector_builder_a2m_input &input,
		index_vector_builder &builder,
		std::string const &chrom_id,
		std::string const &seq_id
	)
	{
		// We don’t assume that the entries appear in particular order, e.g. sorted by chromosome;
		// hence we use only one dispatch group for all operations.
		lb::dispatch_group_async_fn(
			*m_main_group,
			dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),
			[this, bv = builder.destination_bit_vector(), chrom_id, seq_id](){ // Copy everything.
				auto seq_entry(panvc3::msa_index::sequence_entry(seq_id, bv));
				
				{
					std::lock_guard lock{m_msa_index_mutex};
					auto &msa_chr_entry(find_msa_chr_entry(m_msa_index, chrom_id));
					msa_chr_entry.sequence_entries.emplace_back(std::move(seq_entry));
				}
			}
		);
	}
	
	
	void a2m_input_processor::process(input_handler &handler)
	{
		{
			lb::file_ostream msa_index_stream;
			if (!m_msa_index_input_path.empty())
				load_msa_index(m_msa_index_input_path.c_str(), m_msa_index);
		
			lb::open_file_for_writing(m_msa_index_output_path, msa_index_stream, lb::make_writing_open_mode({lb::writing_open_mode::CREATE})); // FIXME: add overwriting conditionally.
			cereal::PortableBinaryOutputArchive msa_archive(msa_index_stream);
		
			lb::log_time(std::cerr) << "Loading the input sequences…\n";
			auto *global_queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
			index_vector_builder_a2m_input a2m_input;
		
			m_main_group.reset(dispatch_group_create());
		
			handler.process_input(m_input_path, [this, &a2m_input](lb::file_handle &handle){
				a2m_input.build(m_index_vector_builder, handle, *this);
			});
			lb::log_time(std::cerr) << "Waiting for compressing the MSA index to finish…\n";
		
			dispatch_group_wait(*m_main_group, DISPATCH_TIME_FOREVER);
			lb::log_time(std::cerr) << "Sorting the index entries…\n";
			std::sort(m_msa_index.chr_entries.begin(), m_msa_index.chr_entries.end());
			for (auto &chr_entry : m_msa_index.chr_entries)
			{
				// Not sure if it’s worth the overhead to parallelise.
				lb::dispatch_group_async_fn(
					*m_main_group,
					dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),
					[&chr_entry](){
						std::sort(chr_entry.sequence_entries.begin(), chr_entry.sequence_entries.end());
					}
				);
			}
		
			dispatch_group_wait(*m_main_group, DISPATCH_TIME_FOREVER);
			lb::log_time(std::cerr) << "Serialising the MSA index…\n";
			msa_archive(m_msa_index);
			msa_index_stream << std::flush;
		}
		
		lb::log_time(std::cerr) << "Done.\n";
		std::exit(EXIT_SUCCESS);
	}
}
