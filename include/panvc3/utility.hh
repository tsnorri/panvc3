/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_UTILITY_HH
#define PANVC3_UTILITY_HH

#include <chrono>
#include <ostream>
#include <seqan3/utility/type_list/type_list.hpp>
#include <tuple>

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
	
	
	template <typename t_value>
	constexpr bool is_power_of_2(t_value const val)
	{
		return 0 < val && !(val & (val - 1));
	}


	// SeqAn 3's program_info_t currently depends on a template parameter, so we need to have a function template for modifying it.
	template <typename t_program_info>
	void append_sam_program_info(
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
		{
			if (pginfo.id.starts_with(id_prefix))
			{
				std::string_view const id(pginfo.id);
				auto const pos(id_prefix.size());
				auto const tail(id.substr(pos));
				std::uint32_t idx_{};
				auto const res(std::from_chars(tail.data(), tail.data() + tail.size(), idx_));
				if (std::errc{} == res.ec)
					idx = std::max(idx, idx_);
			}
		}

		++idx;
		std::string id(id_prefix);
		id += std::to_string(idx);

		// CLI call
		std::string cmd;
		for (int i(0); i < argc; ++i)
		{
			if (0 != i)
				cmd.append(" ");
			cmd.append(argv[i]);
		}

		// Build the record.
		t_program_info rec{
			.id{std::move(id)},
			.command_line_call{std::move(cmd)},
			.version{version}
		};

		rec.name = name;

		if (!dst.empty())
			rec.previous = dst.back().id;

		dst.emplace_back(std::move(rec));
	}


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
}

#endif
