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

class PerfGraphLivePrint : protected ConsoleStreamInterface
{
public:
  PerfGraphLivePrint(PerfGraph & perf_graph, MooseApp & app);

  /// Start printing
  void start();

protected:
  PerfGraph & _perf_graph;

  std::array<PerfGraph::SectionIncrement, MAX_EXECUTION_LIST_SIZE> & _execution_list;

  std::map<PerfID, PerfGraph::SectionInfo> & _id_to_section_info;

  /// This is one beyond the last thing on the stack
  unsigned int _stack_level;

  /// This is one beyond the last thing printed from the stack
  unsigned int _printed_stack_level_end;

  /// The current stack for what the print thread has seen
  std::array<PerfGraph::SectionIncrement, MAX_STACK_SIZE> _print_thread_stack;

  unsigned int _last_execution_list_end;

  bool _printed_name_of_current_section;
};
