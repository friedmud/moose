//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PerfGraph.h"

/**
 * This is effectively a functor that runs on a separate thread and
 * watches the state of the call stack to see if things need to be
 * printed about what the application is doing.
 */
class PerfGraphLivePrint : protected ConsoleStreamInterface
{
public:
  PerfGraphLivePrint(PerfGraph & perf_graph, MooseApp & app);

  /// Start printing
  void start();

protected:
  /**
   * Print the live message
   */
  void printLiveMessage(PerfGraph::SectionIncrement & section_increment);

  /**
   * Print the stats
   */
  void printStats(PerfGraph::SectionIncrement & section_increment_start,
                  PerfGraph::SectionIncrement & section_increment_finish);

  /**
   * Print everything in the stack
   */
  void printStack();

  /**
   * Print everything underneath the current top of the stack
   */
  void printStackUpToLast();

  /**
   * What to do if we're still in the same spot
   */
  void inSamePlace();

  /**
   * What to do if there are new things in the execution list
   */
  void iterateThroughExecutionList();

  /// Number of columns before wrapping
  const unsigned int WRAP_LENGTH = 90;

  /// Reference to the PerfGraph to work with
  PerfGraph & _perf_graph;

  /// Convenience reference to the execution_list within the PerfGraph
  std::array<PerfGraph::SectionIncrement, MAX_EXECUTION_LIST_SIZE> & _execution_list;

  /// This will be true when we should stop printing
  std::future<bool> _done_future;

  /// True when we should stop printing
  const std::atomic<bool> & _destructing;

  /// Whether or not printing is currently turned on
  /// This shadows PerfGraph._live_print_active so that we have consistency
  /// during a single printing
  bool _should_print;

  /// Convenience from PerfGraph
  std::unordered_map<PerfID, PerfGraph::SectionInfo> & _id_to_section_info;

  /// Limit (in seconds) before printing
  Real & _time_limit;

  /// Limit (in MB)
  unsigned int & _mem_limit;

  /// This is one beyond the last thing on the stack
  unsigned int _stack_level;

  /// The current stack for what the print thread has seen
  std::array<PerfGraph::SectionIncrement *, MAX_STACK_SIZE> _print_thread_stack;

  /// The end of the execution list
  /// This is (safely) copied from PerfGraph so that it is consistent for an
  /// entire iteration
  unsigned int _current_execution_list_end;

  /// The actual last entry in the list
  /// This is useful because it is a circular queue - so this is not just
  /// end - 1 (although, often it will be)
  unsigned int _current_execution_list_last;

  /// Where the end of the execution list was during the last call
  /// If this == the current_end... then nothing has happened
  unsigned int _last_execution_list_end;

  /// Which increment was last printed
  PerfGraph::SectionIncrement * _last_printed_increment;

  /// The output count from the console the last time we printed
  unsigned long long int _last_num_printed;

  /// Whether or not printing happened in this iteration
  bool _printed;
};
