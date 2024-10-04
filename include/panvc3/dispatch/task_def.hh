/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_TASK_DEF_HH
#define PANVC3_DISPATCH_TASK_DEF_HH

#include <panvc3/dispatch/queue.hh>
#include <panvc3/dispatch/task_decl.hh>
#include <stdexcept>					// std::logic_error


namespace panvc3::dispatch::detail {
	
	template <typename t_type>
	using add_pointer_to_const_t = std::add_pointer_t <std::add_const_t <t_type>>;
	
	
	template <typename>
	struct is_const_member_function_ {};
	
	template <typename t_return, typename t_class, typename... t_args>
	struct is_const_member_function_ <t_return (t_class::*)(t_args...) const> : std::true_type {};
	
	template <typename t_return, typename t_class, typename... t_args>
	struct is_const_member_function_ <t_return (t_class::*)(t_args...)> : std::false_type {};
	
	template <auto t_fn>
	using is_const_member_function = is_const_member_function_ <decltype(t_fn)>;
	
	template <auto t_fn>
	constexpr static inline auto const is_const_member_function_v{is_const_member_function <t_fn>::value};
	
	
	template <typename t_fn>
	using lambda_target_t = std::conditional_t <
		is_const_member_function_v <&t_fn::operator()>,
		add_pointer_to_const_t <t_fn>,
		std::add_pointer_t <t_fn>
	>;
	
	
	template <typename ... t_args>
	void empty_callable <t_args...>::enqueue_transient_async(queue &qq)
	{
		if constexpr (0 == sizeof...(t_args))
			qq.async(task{empty_callable{}});
		else
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
	}
	
	
	template <typename t_target, typename ... t_args, indirect_member_callable_fn_argt_t <t_target, std::tuple <t_args...>> t_fn>
	void indirect_member_callable <t_target, std::tuple <t_args...>, t_fn>::enqueue_transient_async(queue &qq)
	{
		if constexpr (0 == sizeof...(t_args))
		{
			qq.async(indirect_member_callable <transient_indirect_target_t <t_target>, std::tuple <t_args...>, t_fn>(
				make_transient_indirect_target(target)
			));
		}
		else
		{
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
		}
	}
	
	
	template <typename t_fn, typename ... t_args>
	void lambda_callable <t_fn, t_args...>::enqueue_transient_async(queue &qq)
	{
		if constexpr (0 == sizeof...(t_args))
			qq.async(indirect_member_callable <lambda_target_t <t_fn>, std::tuple <t_args...>, &t_fn::operator()>(&fn));
		else
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
	}
	
	
	template <typename t_fn, typename ... t_args>
	void lambda_ptr_callable <t_fn, t_args...>::enqueue_transient_async(queue &qq)
	{
		if constexpr (0 == sizeof...(t_args))
			qq.async(indirect_member_callable <lambda_target_t <t_fn>, std::tuple <t_args...>, &t_fn::operator()>(fn.get()));
		else
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
	}
	
	
	template <typename ... t_args>
	template <typename t_callable, typename ... t_callable_args>
	task <t_args...>::task(callable_type_tag <t_callable>, t_callable_args && ... args)
	{
		static_assert(sizeof(t_callable) <= BUFFER_SIZE);
		new(m_buffer) t_callable(std::forward <t_callable_args>(args)...);
	}


	// Useful especially for lambda_callable since its exact type can be difficult to pass using callable_type_tag.
	template <typename ... t_args>
	template <typename t_callable>
	task <t_args...>::task(callable_tag, t_callable &&callable)
	{
		static_assert(sizeof(t_callable) <= BUFFER_SIZE);
		new(m_buffer) t_callable(std::forward <t_callable>(callable));
	}
	
	
	template <typename ... t_args>
	auto task <t_args...>::operator=(task &&other) & -> task &
	{
		std::byte buffer[BUFFER_SIZE];
		to_callable(other.m_buffer).move_to(buffer);
		to_callable(m_buffer).move_to(other.m_buffer);
		to_callable(buffer).move_to(m_buffer);
		return *this;
	}
	
	
	template <typename t_target, typename ... t_args, indirect_member_callable_fn_argt_t <t_target, std::tuple <t_args...>> t_fn>
	template <bool>
	void indirect_member_callable <t_target, std::tuple <t_args...>, t_fn>::do_execute(t_args && ... args) requires is_weak_ptr_v <t_target>
	{
		auto target_(target.lock());
		if (target_)
			((*target_).*t_fn)(std::forward <t_args>(args)...);
	}
}

#endif
