/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_TASK_DEF_HH
#define PANVC3_DISPATCH_TASK_DEF_HH

#include <panvc3/dispatch/queue.hh>
#include <panvc3/dispatch/task_decl.hh>
#include <stdexcept>					// std::logic_error


namespace panvc3::dispatch {
	
	template <typename ... t_args>
	void parametrised <t_args...>::empty_callable::enqueue_transient_async(queue &qq)
	{
		if constexpr (0 == sizeof...(t_args))
			qq.async(task{empty_callable{}});
		else
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
	}
	
	
	template <typename ... t_args>
	template <typename t_target, parametrised <t_args...>::indirect_member_callable_fn_t <t_target> t_fn>
	void parametrised <t_args...>::indirect_member_callable <t_target, t_fn>::enqueue_transient_async(queue &qq)
	{
		if constexpr (0 == sizeof...(t_args))
		{
			qq.async(indirect_member_callable <typename transient_indirect_target::type, t_fn>(
				transient_indirect_target::make(target)
			));
		}
		else
		{
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
		}
	}
	
	
	template <typename ... t_args>
	template <typename t_fn>
	void parametrised <t_args...>::lambda_callable <t_fn>::enqueue_transient_async(queue &qq)
	{
		if constexpr (0 == sizeof...(t_args))
			qq.async(indirect_member_callable <t_fn>(&fn));
		else
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
	}
	
	
	template <typename ... t_args>
	template <typename t_fn>
	void parametrised <t_args...>::lambda_ptr_callable <t_fn>::enqueue_transient_async(queue &qq)
	{
		if constexpr (0 == sizeof...(t_args))
			qq.async(indirect_member_callable <t_fn>(fn.get()));
		else
			throw std::logic_error{"enqueue_transient_async() called for callable with parameters"};
	}
	
	
	template <typename ... t_args>
	parametrised <t_args...>::task::task(empty_callable const &)
	{
		static_assert(sizeof(empty_callable) <= BUFFER_SIZE);
		new(m_buffer) empty_callable();
	}
	
	
	template <typename ... t_args>
	template <ClassIndirectCallableTarget <typename parametrised <t_args...>::task, t_args...> t_target>
	parametrised <t_args...>::task::task(t_target &&target)
	{
		typedef indirect_member_callable <std::remove_cvref_t <t_target>> callable_type; // operator() by default.
		static_assert(sizeof(callable_type) <= BUFFER_SIZE);
		new(m_buffer) callable_type(std::forward <t_target>(target));
	}
	
	
	template <typename ... t_args>
	template <
		ClassIndirectTarget <typename parametrised <t_args...>::task, t_args...> t_target,
		parametrised <t_args...>::indirect_member_callable_fn_t <t_target> t_fn,
		typename t_callable
	>
	parametrised <t_args...>::task::task(t_callable &&cc)
	{
		typedef std::remove_cvref_t <t_callable> callable_type;
		static_assert(sizeof(callable_type) <= BUFFER_SIZE);
		new(m_buffer) callable_type(std::forward <t_callable>(cc));
	}
	
	
	template <typename ... t_args>
	template <Callable <t_args...> t_callable>
	parametrised <t_args...>::task::task(t_callable &&callable)
	{
		typedef std::remove_cvref_t <t_callable> callable_type;
		static_assert(sizeof(callable_type) <= BUFFER_SIZE);
		new(m_buffer) callable_type(std::forward <t_callable>(callable));
	}
	
	
	template <typename ... t_args>
	auto parametrised <t_args...>::task::operator=(task &&other) & -> task &
	{
		std::byte buffer[BUFFER_SIZE];
		to_callable(other.m_buffer).move_to(buffer);
		to_callable(m_buffer).move_to(other.m_buffer);
		to_callable(buffer).move_to(m_buffer);
		return *this;
	}
	
	
	template <typename ... t_args>
	template <
		typename t_element,
		parametrised <t_args...>::indirect_member_callable_fn_t <t_element *> t_fn
	>
	template <typename ... t_args_>
	void parametrised <t_args...>::execute_ <std::weak_ptr <t_element>, t_fn>::execute(std::weak_ptr <t_element> &target, t_args_ && ... args)
	{
		auto target_(target.lock());
		if (target_)
			((*target_).*t_fn)(std::forward <t_args_...>(args)...);
	}
}

#endif
