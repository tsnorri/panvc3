/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <filesystem>
#include <iostream>
#include <libbio/dispatch.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/subprocess.hh>
#include <panvc3/msa_index.hh>
#include <string>
#include <string_view>
#include <vector>
#include "cmdline.h"

namespace lb	= libbio;
namespace fs	= std::filesystem;


namespace {
	
	struct chr_entry
	{
		std::string					chr_id;
		std::vector <std::string>	paths;
		
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
	
	
	struct sigchld_handler final : public lb::sigchld_handler
	{
		void child_did_exit_with_nonzero_status(pid_t const pid, int const exit_status, char const *reason) override
		{
			std::cerr << "ERROR: Child process " << pid << " exited with status " << exit_status;
			if (reason)
				std::cerr << " (" << reason << ')';
			std::cerr << '.' << std::endl;
		}
		
		void child_received_signal(pid_t const pid, int const signal_number) override
		{
			std::cerr << "ERROR: Child process " << pid << " received signal " << signal_number << '.' << std::endl;
		}
		
		void finish_handling(bool const did_report_error) override
		{
			if (did_report_error)
				std::exit(EXIT_FAILURE);
		}
	};
	
	
	void install_sigchld_handler()
	{
		static dispatch_once_t once_token;
		dispatch_once(&once_token, ^{
			static sigchld_handler handler;
			lb::install_dispatch_sigchld_handler(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), handler);
		});
	}
	
	
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
	
	
	void read_input_entry(std::size_t const lineno, std::string const &entry, std::string &chr_id, std::string &path)
	{
		std::string_view const &entry_sv(entry);
		std::size_t start(0);
		
		{
			auto const tab_pos(entry_sv.find('\t', start));
			if (std::string::npos == tab_pos)
			{
				std::cerr << "ERROR: Parse error in input on line " << lineno << ".\n";
				std::exit(EXIT_FAILURE);
			}
			
			chr_id = entry_sv.substr(start, tab_pos - start);
			start = 1 + tab_pos;
		}
		
		path = entry_sv.substr(start);
	}


	void load_msa_index(char const *path, panvc3::msa_index &msa_index)
	{
		lb::log_time(std::cerr) << "Loading the input MSA index…\n";
		lb::file_istream stream;
		lb::open_file_for_reading(path, stream);
		cereal::PortableBinaryInputArchive archive(stream);
		archive(msa_index);
	}


	void list_index_contents(char const *path)
	{
		panvc3::msa_index msa_index;
		load_msa_index(path, msa_index);
		
		for (auto const &chr_entry : msa_index.chr_entries)
		{
			std::size_t seq_len{};
			auto const &seq_entries(chr_entry.sequence_entries);
			if (!seq_entries.empty())
				seq_len = seq_entries.front().gap_positions.size();

			std::cout << chr_entry.chr_id << " (" << chr_entry.sequence_entries.size() << " sequences of " << seq_len << " characters)\n";
			for (auto const &seq_entry : seq_entries)
			{
				if (seq_entry.gap_positions.size() != seq_len)
				{
					std::cerr << "ERROR: Sequence " << seq_entry.seq_id << " has " << seq_entry.gap_positions.size() << " characters, expected " << seq_len << '.' << std::endl;
					std::exit(EXIT_FAILURE);
				}

				std::cout << '\t' << seq_entry.seq_id << '\n';
			}
		}
	}


	void query_index(char const *path, char const *chr_id, char const *src_seq_id)
	{
		panvc3::msa_index msa_index;
		load_msa_index(path, msa_index);

		panvc3::msa_index::chr_entry_cmp chr_cmp;
		panvc3::msa_index::sequence_entry_cmp seq_cmp;

		auto const &chr_entries(msa_index.chr_entries);
		auto const chr_rng(std::equal_range(chr_entries.begin(), chr_entries.end(), chr_id, chr_cmp));
		if (chr_rng.first == chr_rng.second)
		{
			std::cerr << "ERROR: No entry for chromosome '" <<  chr_id << "'.\n";
			std::exit(EXIT_FAILURE);
		}

		auto const &chr_entry(*chr_rng.first);
		auto const &seq_entries(chr_entry.sequence_entries);
		auto const seq_rng(std::equal_range(seq_entries.begin(), seq_entries.end(), src_seq_id, seq_cmp));
		if (seq_rng.first == seq_rng.second)
		{
			std::cerr << "ERROR: No entry for sequence '" << src_seq_id << "'.\n";
			std::exit(EXIT_FAILURE);
		}

		auto const &seq_entry(*seq_rng.first);
		auto const aln_limit(seq_entry.gap_positions.size());
		auto const pos_limit(seq_entry.gap_positions_rank0_support(aln_limit));

		while (!std::cin.eof())
		{
			std::cout << "Source co-ordinate? ([0, " << pos_limit << ")) " << std::flush;
			
			std::size_t pos{};
			std::cin >> pos;

			if (! (pos < pos_limit))
				continue;

			std::cout << seq_entry.gap_positions_select0_support(1 + pos) << '\n';
		}
	}


	template <bool t_should_output_seq>
	void build_index(
		std::string const &chr_id,
		std::string const &seq_id,
		lb::file_handle &read_handle,
		std::vector <char> &buffer,
		sdsl::bit_vector &bv,
		std::size_t const wrap_amt
	)
	{
		// Read the input and output non-gapped fasta.
		
		// Prepare for reading.
		buffer.clear();
		struct stat sb;
		read_handle.stat(sb);
		buffer.resize(sb.st_blksize ?: 4096, 0);
		
		bv.clear();
		bv.resize(sb.st_size, 0);
		
		// Output the FASTA header for this sequence.
		if constexpr (t_should_output_seq)
			std::cout << '>' << chr_id << '/' << seq_id << '\n';
		
		// Read and output the sequence.
		std::size_t bytes_read{};
		std::size_t pos{};
		std::size_t non_gap_count{};
		while ((bytes_read = read_handle.read(buffer.size(), buffer.data())))
		{
			for (std::size_t i(0); i < bytes_read; ++i)
			{
				auto const cc(buffer[i]);
				if ('-' == cc)
					bv[pos + i] = 1;
				else
				{
					if constexpr (t_should_output_seq)
					{
						std::cout << cc;
						++non_gap_count;
						
						if (wrap_amt && 0 == non_gap_count % wrap_amt)
							std::cout << '\n';
					}
				}
			}
			
			pos += bytes_read;
		}
		
		if constexpr (t_should_output_seq)
		{
			if (0 == wrap_amt || 0 != non_gap_count % wrap_amt)
				std::cout << '\n';
		}
	}


	struct input_handler
	{
		virtual ~input_handler() {}

		virtual void process_input(
			std::string const &path,
			std::function <void(lb::file_handle &)> const &cb
		) = 0;
	};


	struct file_input_handler final : public input_handler
	{
		virtual void process_input(
			std::string const &path,
			std::function <void(lb::file_handle &)> const &cb
		) override
		{
			lb::file_handle handle(lb::open_file_for_reading(path));
			cb(handle);
		}
	};


	class subprocess_input_handler final : public input_handler
	{
	public:
		typedef lb::subprocess <lb::subprocess_handle_spec::STDOUT> subprocess_type;

	protected:
		std::vector <std::string>	&m_pipe_command;

	public:
		subprocess_input_handler(std::vector <std::string> &pipe_command):
			m_pipe_command(pipe_command)
		{
		}

		virtual void process_input(
			std::string const &path,
			std::function <void(lb::file_handle &)> const &cb
		) override
		{
			m_pipe_command.back() = path;
			auto subprocess(subprocess_type::subprocess_with_arguments(m_pipe_command));
			cb(subprocess.stdout_handle());
		}
	};
	
	
	class input_processor
	{
	protected:
		std::string					m_input_path;
		std::string					m_msa_index_input_path;
		std::string					m_msa_index_output_path;
		std::vector <std::string>	m_pipe_command;
		std::size_t					m_fasta_line_width{};
		bool						m_should_output_fasta{};
		
	public:
		input_processor(
			char const *input_path,
			char const *msa_index_input_path,
			char const *msa_index_output_path,
			char const *pipe_input_command,
			bool const should_output_fasta,
			std::size_t const fasta_line_width
		):
			m_input_path(input_path),
			m_msa_index_input_path(msa_index_input_path ?: ""),
			m_msa_index_output_path(msa_index_output_path),
			m_pipe_command(lb::parse_command_arguments(pipe_input_command)),
			m_fasta_line_width(fasta_line_width),
			m_should_output_fasta(should_output_fasta)
		{
		}
		

		void process(input_handler &handler)
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
				std::string path;
				std::size_t lineno(1);
				
				while (std::getline(path_stream, entry_buffer))
				{
					read_input_entry(lineno, entry_buffer, chr_id, path);
					auto &entry(find_chr_entry(chr_entries, chr_id));
					entry.paths.emplace_back(path);
					++lineno;
				}
			}
			
			// Prepare the MSA index.
			msa_index.chr_entries.reserve(msa_index.chr_entries.size() + chr_entries.size());
			
			// Handle the inputs and compress and sort in background.
			std::vector <char> seq_buffer;
			lb::dispatch_ptr <dispatch_group_t> main_group(dispatch_group_create());
			auto global_queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
	
			for (auto const &entry : chr_entries)
			{
				lb::log_time(std::cerr) << "Handling sequences for chromosome " << entry.chr_id << "…\n";
				
				lb::dispatch_ptr <dispatch_group_t> chr_group(dispatch_group_create());
				auto &msa_chr_entry(find_msa_chr_entry(msa_index, entry.chr_id));
				msa_chr_entry.sequence_entries.reserve(entry.paths.size());
				
				for (auto const &path_str : entry.paths)
				{
					handler.process_input(
						path_str,
						[
							this,
							global_queue,
							&seq_buffer,
							&entry,
							&chr_group,
							&msa_chr_entry,
							&path_str
						](lb::file_handle &handle){
							sdsl::bit_vector bv;
							fs::path const path(path_str);
							auto const fname(path.filename());
							lb::log_time(std::cerr) << "Processing " << fname << "…\n";

							if (m_should_output_fasta)
								build_index <true>(entry.chr_id, fname, handle, seq_buffer, bv, m_fasta_line_width);
							else
								build_index <false>(entry.chr_id, fname, handle, seq_buffer, bv, m_fasta_line_width);
							
							auto &seq_entry(msa_chr_entry.sequence_entries.emplace_back());
							lb::dispatch_group_async_fn(*chr_group, global_queue, [fname, bv = std::move(bv), &seq_entry](){
								seq_entry = panvc3::msa_index::sequence_entry(std::move(fname), bv);
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


		void operator()()
		{
			if (m_pipe_command.empty())
			{
				file_input_handler handler;
				process(handler);
			}
			else
			{
				subprocess_input_handler handler(m_pipe_command);
				process(handler);
			}
		}
	};
	
	
	extern void process(gengetopt_args_info &args_info)
	{
		if (args_info.list_contents_given)
		{
			list_index_contents(args_info.msa_index_input_arg);
			std::exit(EXIT_SUCCESS);
		}
		else if (args_info.query_given)
		{
			query_index(args_info.msa_index_input_arg, args_info.chr_id_arg, args_info.src_seq_id_arg);
			std::exit(EXIT_SUCCESS);
		}
		else if (args_info.build_index_given)
		{
			static input_processor processor(
				args_info.sequence_inputs_arg,
				args_info.msa_index_input_arg,
				args_info.msa_index_output_arg,
				args_info.pipe_input_arg,
				args_info.output_fasta_flag,
				args_info.fasta_line_width_arg
			);
			
			auto caller(lb::dispatch(processor));
			caller.async <>(dispatch_get_main_queue());
		}
		else
		{
			std::cerr << "ERROR: Unknown mode.\n";
			std::exit(EXIT_FAILURE);
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	
	if (args_info.fasta_line_width_arg < 0)
	{
		std::cerr << "ERROR: FASTA line width must be non-negative.\n";
		std::exit(EXIT_FAILURE);
	}
	
	install_sigchld_handler();
	process(args_info);
	
	dispatch_main();

	// Not reached.
	return EXIT_SUCCESS;
}
