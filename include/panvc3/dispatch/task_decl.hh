/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_TASK_DECL_HH
#define PANVC3_DISPATCH_TASK_DECL_HH

#include <cstddef>					// std::byte
#include <memory>					// std::weak_ptr
#include <panvc3/dispatch/fwd.hh>
#include <utility>					// std::forward


namespace panvc3::dispatch::detail {
	
	template <typename t_ptr>
	struct ptr_value_type { typedef typename t_ptr::element_type type; }; // For smart pointers.
	
	template <typename t_type>
	struct ptr_value_type <t_type *> { typedef t_type type; };
	
	template <typename t_ptr>
	using ptr_value_type_t = typename ptr_value_type <t_ptr>::type;
	
	
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
}


namespace panvc3::dispatch {
	
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
	
	
	template <typename ... t_args>
	struct parametrised
	{
		class task; // Fwd.
		constexpr static inline std::size_t const TASK_BUFFER_SIZE{4 * sizeof(void *)}; // Enough for “most” cases; holds vptr + shared_ptr + void *.
		
		
		// Not in namespace detail b.c. depends on t_targs and not repeating the struct template seemed easier.
		typedef void(*indirect_callable_fn_type)(t_args...); // For non-member functions.
		
		template <typename t_target>
		struct indirect_member_callable_fn
		{
			typedef void(detail::ptr_value_type_t <t_target>::*type)(t_args...); // For member functions.
		};
		
		template <typename t_target>
		using indirect_member_callable_fn_t = typename indirect_member_callable_fn <t_target>::type;
		
		
		template <typename t_target, bool = std::is_invocable_v <t_target, t_args...>> // Check for operator().
		struct indirect_member_callable_default_fn
		{
			constexpr static inline indirect_member_callable_fn_t <t_target> const value{nullptr}; // No default target.
		};
		
		template <typename t_target>
		struct indirect_member_callable_default_fn <t_target, true>
		{
			constexpr static inline indirect_member_callable_fn_t <t_target> const value{&t_target::operator()}; // Call operator().
		};
		
		template <typename t_target>
		constexpr static inline indirect_member_callable_fn_t <t_target> indirect_member_callable_default_fn_v = indirect_member_callable_default_fn <t_target>::value;
		
		
		// Import here.
		typedef detail::callable <t_args...>	callable;
		
		struct empty_callable final : public callable
		{
			void execute(t_args && ...) override {}
			void move_to(std::byte *buffer) override { new(buffer) empty_callable{}; }
			inline void enqueue_transient_async(queue &qq) override;
		};
		
		
		template <typename t_fn>
		struct lambda_callable final : public callable
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
		
		
		template <typename t_fn>
		struct lambda_ptr_callable final : public callable
		{
			std::unique_ptr <t_fn>	fn;
			
			lambda_ptr_callable(std::unique_ptr <t_fn> &&fn_):
				fn(std::move(fn_))
			{
			}
			
			lambda_ptr_callable(t_fn &&fn_) requires std::is_rvalue_reference_v <t_fn>:
				lambda_ptr_callable(std::make_unique(std::move(fn_)))
			{
			}
			
			void execute(t_args && ... args) override { (*fn)(std::forward <t_args>(args)...); }
			void move_to(std::byte *buffer) override { new(buffer) lambda_ptr_callable{std::move(fn)}; }
			inline void enqueue_transient_async(queue &qq) override;
		};
		
		
		template <typename t_fn>
		static lambda_callable <t_fn> make_lambda_callable(t_fn &&fn)
		requires (sizeof(lambda_callable <t_fn>) <= TASK_BUFFER_SIZE)
		{
			return lambda_callable{std::forward <t_fn>(fn)};
		}
		
		template <typename t_fn>
		static lambda_ptr_callable <t_fn> make_lambda_callable(t_fn &&fn)
		requires (TASK_BUFFER_SIZE < sizeof(lambda_callable <t_fn>))
		{
			return lambda_ptr_callable{std::forward <t_fn>(fn)};
		}
		
		
		// FIXME: add a subclass for non-member functions.
		
		
		template <typename t_target, indirect_member_callable_fn_t <t_target> t_fn>
		struct execute_
		{
			template <typename ... t_args_>
			static void execute(t_target &target, t_args_ && ... args) { ((*target).*t_fn)(std::forward <t_args_...>(args)...); }
		};
		
		template <typename t_element, indirect_member_callable_fn_t <t_element *> t_fn>
		struct execute_ <std::weak_ptr <t_element>, t_fn>
		{
			template <typename ... t_args_>
			static inline void execute(std::weak_ptr <t_element> &target, t_args_ && ... args);
		};
		
		
		// Small b.c. the member function pointer is passed as a template parameter.
		template <typename t_target, indirect_member_callable_fn_t <t_target> t_fn = indirect_member_callable_default_fn_v <t_target>>
		struct indirect_member_callable final : public callable
		{
			typedef detail::transient_indirect_target <t_target>	transient_indirect_target;
			
			t_target	target;
			
			template <typename t_target_>
			indirect_member_callable(t_target_ &&target_):
				target(std::forward <t_target_>(target_))
			{
			}
			
			void execute(t_args && ... args) override { execute_ <t_target, t_fn>::execute(target, std::forward <t_args...>(args)...); }
			void move_to(std::byte *buffer) override { new(buffer) indirect_member_callable{std::move(target)}; }
			inline void enqueue_transient_async(queue &qq) override;
		};
		
		
		class alignas(alignof(std::max_align_t)) task // Take into account the alignment of every possible t_fn and t_target.
		{
		private:
			constexpr static inline std::size_t const BUFFER_SIZE{TASK_BUFFER_SIZE};
			static callable &to_callable(std::byte *buffer) { return *reinterpret_cast <callable *>(buffer); }
			
		private:
			std::byte m_buffer[BUFFER_SIZE];
			
		private:
			template <
				ClassIndirectTarget <task, t_args...> t_target,
				indirect_member_callable_fn_t <t_target> t_fn,
				typename t_callable
			>
			explicit inline task(t_callable &&callable);
			callable &get_callable() { return to_callable(m_buffer); }
			
		public:
			virtual ~task() { get_callable().~callable(); }
			
			// Construct from an empty_callable.
			/* implicit */ inline task(empty_callable const &);
			
			// Default constructor, produces an empty_callable.
			task(): task(empty_callable{}) {}
			
			// Construct from a target that produces an indirect_member_callable. (I.e. pointer, std::unique_ptr, std::shared_ptr, std::weak_ptr.)
			template <ClassIndirectCallableTarget <task, t_args...> t_target>
			/* implicit */ inline task(t_target &&target);
			
			// Construct from a callable.
			template <Callable <t_args...> t_callable>
			/* implicit */ inline task(t_callable &&cc);
			
			// Construct from a target that produces an indirect_member_callable. (I.e. pointer, std::unique_ptr, std::shared_ptr, std::weak_ptr.)
			// The size check is in the constructor.
			template <ClassIndirectTarget <task, t_args...> t_target, indirect_member_callable_fn_t <t_target> t_fn>
			static task from_member_fn(t_target &&target) { return task{indirect_member_callable <t_target, t_fn>{std::forward <t_target>(target)}}; }
			
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
	};
	
	
	typedef parametrised <>::task task;
}

#endif
