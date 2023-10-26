/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_TASK_DECL_HH
#define PANVC3_DISPATCH_TASK_DECL_HH

#include <cstddef>					// std::byte
#include <libbio/tuple/access.hh>	// libbio::tuples::forward_to
#include <memory>					// std::weak_ptr
#include <panvc3/dispatch/fwd.hh>
#include <utility>					// std::forward


namespace panvc3::dispatch::detail {
	
	template <typename...> struct callable; // Fwd.
}


namespace panvc3::dispatch {
	
	constexpr static inline std::size_t const TASK_BUFFER_SIZE{4 * sizeof(void *)}; // Enough for “most” cases; holds vptr + shared_ptr + void *.
	
	// We would like to distinguish the types whose base is callable but not erase their type.
	template <typename t_type, typename ... t_args>
	concept Callable = std::is_base_of_v <detail::callable <t_args...>, std::remove_cvref_t <t_type>>;
	
	// Types that are:
	// – Not task
	// – Not Callable
	template <typename t_target, typename t_task, typename ... t_args>
	concept ClassIndirectTarget =
		// FIXME: Check that t_target is class? Should get the underlying type from different smart pointers.
 		(!std::is_same_v <t_task, std::remove_cvref_t <t_target>>) &&
		(!Callable <t_target, t_args...>);
	
	// Do not inherit callable (and hence not Callable) b.c. the target is accessed via a pointer;
	// no need to move etc.
	// Also have the function call operator via pointer.
	template <typename t_target, typename t_task, typename ... t_args>
	concept ClassIndirectCallableTarget =
		ClassIndirectTarget <t_target, t_task, t_args...> &&
		requires (t_target target, t_args... args)
		{
			(*target)(args...);
		};
}


namespace panvc3::dispatch::detail {
	
	template <typename t_ptr>
	struct ptr_value_type { typedef typename t_ptr::element_type type; }; // For smart pointers.
	
	template <typename t_type>
	struct ptr_value_type <t_type *> { typedef t_type type; };
	
	template <typename t_ptr>
	using ptr_value_type_t = typename ptr_value_type <t_ptr>::type;
	
	
	// Type trait for arguments passed in a tuple.
	template <typename, typename>
	struct is_invocable {};
	
	template <typename t_target, typename ... t_args>
	struct is_invocable <t_target, std::tuple <t_args...>> : public std::is_invocable <t_target, t_args...> {};
	
	template <typename t_target, typename t_args>
	constexpr static inline auto const is_invocable_v{is_invocable <t_target, t_args>::value};
	
	
	// Type trait for std::weak_ptr
	template <typename>
	struct is_weak_ptr : public std::false_type {};
	
	template <typename t_type>
	struct is_weak_ptr <std::weak_ptr <t_type>> : public std::true_type {};
	
	template <typename t_type>
	constexpr static inline bool const is_weak_ptr_v{is_weak_ptr <t_type>::value};
	
	
	// A type that can call different types of functions, including member functions (in classes with virtual base classes).
	// Here b.c. we would like to reference it in a concept.
	template <typename ... t_args>
	struct callable
	{
		virtual ~callable() {}
		virtual void execute(t_args && ...) = 0;
		virtual void move_to(std::byte *buffer) = 0;
		virtual void enqueue_transient_async(queue &qq) = 0; // Make a trasient version of the callable and execute it asynchronously in the given queue; used by event sources.
		
		template <typename ... t_args_>
		void operator()(t_args_ && ... args) { execute(std::forward <t_args_...>(args)...); }
	};
	
	
	template <typename ... t_args>
	using indirect_callable_fn_type = void (*)(t_args...); // For non-member functions.
	
	
	template <typename t_target, typename... t_args>
	struct indirect_member_callable_fn
	{
		typedef ptr_value_type_t <t_target> value_type;
		typedef std::conditional_t <
			std::is_const_v <value_type>,
			void(value_type::*)(t_args...) const,
			void(value_type::*)(t_args...)
		> type; // For member functions.
	};
	
	template <typename t_target, typename ... t_args>
	using indirect_member_callable_fn_t = typename indirect_member_callable_fn <t_target, t_args...>::type;
	
	
	// Take the parameters in a tuple.
	template <typename, typename>
	struct indirect_member_callable_fn_argt {};
	
	template <typename t_target, typename... t_args>
	struct indirect_member_callable_fn_argt <t_target, std::tuple <t_args...>>
	{
		typedef indirect_member_callable_fn_t <t_target, t_args...> type;
	};
	
	template <typename t_target, typename t_arg_tuple>
	using indirect_member_callable_fn_argt_t = typename indirect_member_callable_fn_argt <t_target, t_arg_tuple>::type;
	
	
	template <typename t_target, typename t_arg_tuple, bool t_is_invocable = is_invocable_v <t_target, t_arg_tuple>> // Check for operator().
	struct indirect_member_callable_default_fn
	{
		template <typename ... t_args_>
		using fn_t = indirect_member_callable_fn_t <t_target, t_args_...>;
		
		typedef libbio::tuples::forward_to <fn_t>::template parameters_of_t <t_arg_tuple> fn_type;
		
		constexpr static inline fn_type const value{nullptr}; // No default target.
	};
	
	template <typename t_target, typename t_arg_tuple>
	struct indirect_member_callable_default_fn <t_target, t_arg_tuple, true>
	{
		template <typename ... t_args_>
		using fn_t = indirect_member_callable_fn_t <t_target, t_args_...>;
		
		typedef libbio::tuples::forward_to <fn_t>::template parameters_of_t <t_arg_tuple> fn_type;
		
		constexpr static inline fn_type const value{&t_target::operator()}; // Call operator().
	};
	
	template <typename t_target, typename t_arg_tuple>
	constexpr static inline auto const indirect_member_callable_default_fn_argt_v{
		indirect_member_callable_default_fn <t_target, t_arg_tuple>::value
	};
	
	template <typename t_target, typename ... t_args>
	constexpr static inline auto const indirect_member_callable_default_fn_v{
		indirect_member_callable_default_fn_argt_v <t_target, std::tuple <t_args...>>
	};
	
	
	template <typename t_target>
	struct transient_indirect_target
	{
		typedef t_target	type;
		static type &make(type &tt) { return tt; }
	};
	
	template <typename t_element>
	struct transient_indirect_target <std::unique_ptr <t_element>>
	{
		typedef t_element	*type;
		static type make(std::unique_ptr <t_element> &ptr) { return ptr.get(); }
	};
	
	template <typename t_type>
	using transient_indirect_target_t = typename transient_indirect_target <t_type>::type;
	
	template <typename t_target>
	auto make_transient_indirect_target(t_target &&target) { return transient_indirect_target <t_target>::make(target); }
	
	
	template <typename ... t_args>
	struct empty_callable final : public callable <t_args...>
	{
		void execute(t_args && ...) override {}
		void move_to(std::byte *buffer) override { new(buffer) empty_callable{}; }
		inline void enqueue_transient_async(queue &qq) override;
	};
	
	
	template <typename t_fn, typename ... t_args>
	struct lambda_callable final : public callable <t_args...>
	{
		t_fn	fn;
		
		lambda_callable(t_fn &&fn_): //requires std::is_rvalue_reference_v <t_fn>:
			fn(std::move(fn_))
		{
		}
		
		void execute(t_args && ... args) override { fn(std::forward <t_args>(args)...); }
		void move_to(std::byte *buffer) override { new(buffer) lambda_callable{std::move(fn)}; }
		inline void enqueue_transient_async(queue &qq) override;
	};
	
	
	template <typename t_fn, typename ... t_args>
	struct lambda_ptr_callable final : public callable <t_args...>
	{
		std::unique_ptr <t_fn>	fn;
		
		lambda_ptr_callable(std::unique_ptr <t_fn> &&fn_):
			fn(std::move(fn_))
		{
		}
		
		lambda_ptr_callable(t_fn &&fn_): //requires std::is_rvalue_reference_v <t_fn>:
			lambda_ptr_callable(std::make_unique <t_fn>(std::move(fn_)))
		{
		}
		
		void execute(t_args && ... args) override { (*fn)(std::forward <t_args>(args)...); }
		void move_to(std::byte *buffer) override { new(buffer) lambda_ptr_callable{std::move(fn)}; }
		inline void enqueue_transient_async(queue &qq) override;
	};
	
	
	template <typename t_fn, typename... t_args>
	static lambda_callable <t_fn, t_args...> make_lambda_callable(t_fn &&fn)
	requires (sizeof(lambda_callable <t_fn, t_args...>) <= TASK_BUFFER_SIZE)
	{
		return lambda_callable <t_fn, t_args...>{std::forward <t_fn>(fn)};
	}
	
	template <typename t_fn, typename... t_args>
	static lambda_ptr_callable <t_fn, t_args...> make_lambda_callable(t_fn &&fn)
	requires (TASK_BUFFER_SIZE < sizeof(lambda_callable <t_fn, t_args...>))
	{
		return lambda_ptr_callable <t_fn, t_args...>{std::forward <t_fn>(fn)};
	}
	
	
	// FIXME: add a subclass for non-member functions.
	
	
	// Small b.c. the member function pointer is passed as a template parameter.
	template <typename t_target, typename t_arg_tuple, indirect_member_callable_fn_argt_t <t_target, t_arg_tuple> t_fn>
	struct indirect_member_callable {};
	
	template <typename t_target, typename ... t_args, indirect_member_callable_fn_argt_t <t_target, std::tuple <t_args...>> t_fn>
	struct indirect_member_callable <t_target, std::tuple <t_args...>, t_fn> : public callable <t_args...>
	{
		t_target	target;
		
		template <typename t_target_>
		indirect_member_callable(t_target_ &&target_):
			target(std::forward <t_target_>(target_))
		{
		}
		
	private:
		template <bool = true>
		inline void do_execute(t_args && ... args) requires is_weak_ptr_v <t_target>;
		
		template <bool = true>
		void do_execute(t_args && ... args) requires (!is_weak_ptr_v <t_target>) { ((*target).*t_fn)(std::forward <t_args>(args)...); }
	
	public:
		void execute(t_args && ... args) override { do_execute(std::forward <t_args>(args)...); }
		
		void move_to(std::byte *buffer) override { new(buffer) indirect_member_callable{std::move(target)}; }
		inline void enqueue_transient_async(queue &qq) override;
	};
	
	
	template <typename... t_args>
	class alignas(alignof(std::max_align_t)) task // Take into account the alignment of every possible t_fn and t_target.
	{
	public:
		typedef callable <t_args...>		callable_type;
		typedef empty_callable <t_args...>	empty_callable_type;
		
		template <typename t_target>
		constexpr static inline auto const indirect_member_callable_default_fn_v{
			detail::indirect_member_callable_default_fn_v <t_target, t_args...>
		};
		
		template <typename t_target, indirect_member_callable_fn_t <t_target, t_args...> t_fn = indirect_member_callable_default_fn_v <t_target>>
		using indirect_member_callable_t = indirect_member_callable <t_target, std::tuple <t_args...>, t_fn>;
		
	private:
		constexpr static inline std::size_t const BUFFER_SIZE{TASK_BUFFER_SIZE};
		static callable_type &to_callable(std::byte *buffer) { return *reinterpret_cast <callable_type *>(buffer); }
		
		template <typename t_fn>
		static auto make_lambda_callable(t_fn &&fn) { return detail::make_lambda_callable <t_fn, t_args...>(std::forward <t_fn>(fn)); }
		
	private:
		std::byte m_buffer[BUFFER_SIZE];
		
	private:
		template <typename t_type>
		struct private_tag {};
		
		template <typename t_callable, typename ... t_args_>
		explicit inline task(private_tag <t_callable>, t_args_ && ... args);
		
		callable_type &get_callable() { return to_callable(m_buffer); }
		
	public:
		virtual ~task() { get_callable().~callable(); }
		
		// Default constructor, produces an empty_callable.
		task(): task{private_tag <empty_callable_type>{}} {}
		
		// Construct from an empty_callable.
		/* implicit */ inline task(empty_callable_type const &): task{private_tag <empty_callable_type>{}} {}
		
		// Construct from a target that produces an indirect_member_callable. (I.e. pointer, std::unique_ptr, std::shared_ptr, std::weak_ptr.)
		// (See note above about task <t_args...>.)
		template <ClassIndirectCallableTarget <task <t_args...>, t_args...> t_target>
		/* implicit */ inline task(t_target &&target):
			task{private_tag <indirect_member_callable_t <std::remove_reference_t <t_target>>>{}, std::forward <t_target>(target)} {} // operator() by default
		
		// Construct from a callable.
		template <Callable <t_args...> t_callable>
		/* implicit */ inline task(t_callable &&callable):
			task{private_tag <std::remove_reference_t <t_callable>>{}, std::forward <t_callable>(callable)} {}
		
		// Construct from a target that produces an indirect_member_callable. (I.e. pointer, std::unique_ptr, std::shared_ptr, std::weak_ptr.)
		// The size check is in the constructor.
		// (See note above about task <t_args...>.)
		template <ClassIndirectTarget <task <t_args...>, t_args...> t_target, indirect_member_callable_fn_t <t_target, t_args...> t_fn>
		static task from_member_fn(t_target &&target) { return task{private_tag <indirect_member_callable_t <t_target, t_fn>>{}, std::forward <t_target>(target)}; }
		
		template <typename t_fn>
		static task from_lambda(t_fn &&fn) { return make_lambda_callable(std::forward <t_fn>(fn)); }
		
		// Move constructor and assignment.
		task(task &&other) { other.get_callable().move_to(m_buffer); }
		inline task &operator=(task &&other) &;
		
		// Likely no need to copy.
		task(task const &) = delete;
		task &operator=(task const &) = delete;
		
		template <typename ... t_args_>
		void execute(t_args_ && ... args) { get_callable().execute(std::forward <t_args_...>(args)...); }
		
		template <typename ... t_args_>
		void operator()(t_args_ && ... args) { execute(std::forward <t_args_...>(args)...); }
		
		void enqueue_transient_async(queue &qq) requires (0 == sizeof...(t_args)) { get_callable().enqueue_transient_async(qq); }
	};
}


namespace panvc3::dispatch {
	
	template <typename ... t_args>
	using task_t = detail::task <t_args...>;
	
	typedef detail::task <> task;
}

#endif
