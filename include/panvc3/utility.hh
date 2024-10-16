/*
 * Copyright (c) 2022-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_UTILITY_HH
#define PANVC3_UTILITY_HH

#include <chrono>
#include <condition_variable>
#include <libbio/dispatch.hh>
#include <libbio/sam/header.hh>
#include <mutex>
#include <ostream>
#include <tuple>

#if defined(PANVC3_USE_SEQAN3) && PANVC3_USE_SEQAN3
#include <seqan3/utility/type_list/type_list.hpp>
#endif


namespace panvc3 {

	// FIXME: Move to libbio?
	template <typename t_type, typename t_proj, typename t_cmp = std::less <>>
	struct cmp_proj
	{
	private:
		template <typename, typename = void> struct is_transparent_t : public std::false_type {};

		template <typename t_type_>
		struct is_transparent_t <t_type_, std::void_t <typename t_type_::is_transparent>> : public t_type_::is_transparent {};

	public:
		typedef is_transparent_t <t_cmp>	is_transparent;


		constexpr bool operator()(t_type const &lhs, t_type const &rhs) const
		{
			t_cmp cmp;
			t_proj proj;
			return cmp(proj(lhs), proj(rhs));
		}

		template <typename t_other>
		constexpr bool operator()(t_type const &lhs, t_other const &rhs) const
		{
			t_cmp cmp;
			t_proj proj;
			return cmp(proj(lhs), rhs);
		}

		template <typename t_other>
		constexpr bool operator()(t_other const &lhs, t_type const &rhs) const
		{
			t_cmp cmp;
			t_proj proj;
			return cmp(lhs, proj(rhs));
		}
	};
	

	// FIXME: check that t_alphabet is one of SeqAn3’s alphabets.
	template <typename t_alphabet>
	constexpr inline t_alphabet max_letter()
	{
		t_alphabet retval;
		retval.assign_rank(t_alphabet::alphabet_size - 1);
		return retval;
	}
	
	
	std::uint32_t check_sam_program_id_index(std::string_view const id_prefix, std::string_view const id);
	
	std::string command_line_call(
		int const argc,
		char const * const * const argv
	);
	
	void append_sam_program_info(
		std::string_view const id_prefix,
		std::string_view const name,
		std::string call,
		char const * const version,
		libbio::sam::program_entry_vector &dst
	);
	
	void append_sam_program_info(
		std::string_view const id_prefix,
		std::string_view const name,
		int const argc,
		char const * const * const argv,
		char const * const version,
		libbio::sam::program_entry_vector &dst
	);
	
	
	void prepare_thread_pool_with_args(libbio::dispatch::thread_pool &thread_pool, long threads); // Type from gengetopt
	
	
	template <typename t_duration>
	std::ostream &log_duration(std::ostream &os, t_duration dur)
	{
		namespace chrono = std::chrono;

		typedef chrono::duration <std::uint64_t, std::ratio <3600 * 24>> day_type;
		auto const dd{chrono::duration_cast <day_type>(dur)};
		auto const hh{chrono::duration_cast <chrono::hours>(dur -= dd)};
		auto const mm{chrono::duration_cast <chrono::minutes>(dur -= hh)};
		auto const ss{chrono::duration_cast <chrono::seconds>(dur -= mm)};
		auto const ms{chrono::duration_cast <chrono::milliseconds>(dur -= ss)};

		bool should_print{false};
		auto const dd_(dd.count());
		if (dd_)
		{
			os << dd_ << " d, ";
			should_print = true;
		}

		auto const hh_(hh.count());
		if (should_print || hh_)
		{
			os << hh_ << " h, ";
			should_print = true;
		}

		auto const mm_(mm.count());
		if (should_print || mm_)
		{
			os << mm_ << " m, ";
			should_print = true;
		}

		auto const ss_(ss.count());
		if (should_print || ss_)
			os << ss_ << " s, ";

		os << ms.count() << " ms";

		return os;
	}
	
	
	struct timer
	{
	private:
		mutable std::mutex				m_mutex{};
		mutable std::condition_variable	m_cv{};
		bool							m_should_stop{};

	public:
		inline void stop();

		template <typename t_rep, typename t_period>
		bool wait_for(std::chrono::duration <t_rep, t_period> const duration) const;
	};


	void timer::stop()
	{
		std::unique_lock lock(m_mutex);
		m_should_stop = true;
		m_cv.notify_all();
	}


	template <typename t_rep, typename t_period>
	bool timer::wait_for(std::chrono::duration <t_rep, t_period> const duration) const
	{
		std::unique_lock lock(m_mutex);
		return !m_cv.wait_for(lock, duration, [this](){ return m_should_stop; });
	}
}

#if defined(PANVC3_USE_SEQAN3) && PANVC3_USE_SEQAN3
namespace panvc3 {
	
	// I’m not sure how to get the types from seqan3::type_list, so here is a helper for that.
	// (Also I don’t know if there is syntax for generalising this for any template that takes a parameter pack.)
	
	template <typename t_type>
	struct type_list_to_tuple {};
	
	template <typename... t_types>
	struct type_list_to_tuple <seqan3::type_list <t_types...>>
	{
		typedef std::tuple <t_types...> type;
	};
	
	template <typename t_type>
	using type_list_to_tuple_t = typename type_list_to_tuple <t_type>::type;
	
	
	// SeqAn 3's program_info_t currently depends on a template parameter, so we need to have a function template for modifying it.
	template <typename t_program_info>
	void append_sam_program_info_seqan3(
		std::string_view const id_prefix,
		std::string_view const name,
		int const argc,
		char const * const * const argv,
		char const * const version,
		std::vector <t_program_info> &dst
	)
	{
		// Identifier
		std::uint32_t idx{};
		for (auto const &pginfo : dst)
			idx = std::max(idx, check_program_id_index(id_prefix, pginfo.id));

		++idx;
		std::string id(id_prefix);
		id += std::to_string(idx);

		// Build the record.
		t_program_info rec{
			.id{std::move(id)},
			.command_line_call{command_line_call(argc, argv)},
			.version{version}
		};

		rec.name = name;

		if (!dst.empty())
			rec.previous = dst.back().id;

		dst.emplace_back(std::move(rec));
	}
}
#endif

#endif
