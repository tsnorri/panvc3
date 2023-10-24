/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_HH
#define PANVC3_DISPATCH_HH

// Ownership Semantics
// ===================
// – The library user owns the queues, the groups, the thread pool and the event manager.
// – The event manager owns the event sources.
// – Task ownership depends on the task type.
//
// Other things that are noteworthy
// ================================
// – A barrier is guaranteed to stop a queue cleanly if it calls thread_pool::stop() synchronously on the thread in which its task is executed.

#include <panvc3/dispatch/barrier.hh>
#include <panvc3/dispatch/group.hh>
#include <panvc3/dispatch/queue.hh>
#include <panvc3/dispatch/task_decl.hh>
#include <panvc3/dispatch/task_def.hh>

#endif
