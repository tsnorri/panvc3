/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_TASK_DEF_HH
#define PANVC3_DISPATCH_TASK_DEF_HH

#include <panvc3/dispatch/queue.hh>
#include <panvc3/dispatch/task_decl.hh>
#include <stdexcept>					// std::logic_error


namespace panvc3::dispatch::detail {
	
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
			qq.async(indirect_member_callable <t_fn, std::tuple <t_args...>, &t_fn::operator()>(&fn));
		else
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
	}
	
	
	template <typename t_fn, typename ... t_args>
	void lambda_ptr_callable <t_fn, t_args...>::enqueue_transient_async(queue &qq)
	{
		if constexpr (0 == sizeof...(t_args))
			qq.async(indirect_member_callable <t_fn, std::tuple <t_args...>, &t_fn::operator()>(fn.get()));
		else
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
	}
	
	
	template <typename ... t_args>
	task <t_args...>::task(empty_callable_type const &)
	{
		static_assert(sizeof(empty_callable_type) <= BUFFER_SIZE);
		new(m_buffer) empty_callable();
	}
	
	
	template <typename ... t_args>
	template <ClassIndirectCallableTarget <task <t_args...>, t_args...> t_target>
	task <t_args...>::task(t_target &&target)
	{
		typedef indirect_member_callable_t <std::remove_cvref_t <t_target>> callable_type; // operator() by default.
		static_assert(sizeof(callable_type) <= BUFFER_SIZE);
		new(m_buffer) callable_type(std::forward <t_target>(target));
	}
	
	
	template <typename ... t_args>
	template <
		ClassIndirectTarget <task <t_args...>, t_args...> t_target,
		indirect_member_callable_fn_t <t_target, t_args...> t_fn,
		typename t_callable
	>
	task <t_args...>::task(t_callable &&cc)
	{
		typedef std::remove_cvref_t <t_callable> callable_type;
		static_assert(sizeof(callable_type) <= BUFFER_SIZE);
		new(m_buffer) callable_type(std::forward <t_callable>(cc));
	}
	
	
	template <typename ... t_args>
	template <Callable <t_args...> t_callable>
	task <t_args...>::task(t_callable &&callable)
	{
		typedef std::remove_cvref_t <t_callable> callable_type;
		static_assert(sizeof(callable_type) <= BUFFER_SIZE);
		new(m_buffer) callable_type(std::forward <t_callable>(callable));
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
