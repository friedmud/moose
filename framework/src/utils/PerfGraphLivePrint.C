//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PerfGraphLivePrint.h"

PerfGraphLivePrint::PerfGraphLivePrint(PerfGraph & perf_graph, MooseApp & app) : ConsoleStreamInterface(app), _perf_graph(perf_graph),

                                                                 _execution_list(perf_graph._execution_list),
                                                                 _id_to_section_info(perf_graph._id_to_section_info),
                                                                 _stack_level(0),
                                                                 _printed_stack_level_end(0),
                                                                 _last_execution_list_end(0),
                                                                 _printed_name_of_current_section(false)
{
}

void
PerfGraphLivePrint::start()
{
  while(true)
  {
    std::this_thread::sleep_for(std::chrono::seconds(1));

    // Where we were last time in the execution list
    // this should point at the first _new_ thing... or at end
//      auto current_execution_list_begin = this->_execution_list_begin.load(std::memory_order_relaxed);

    // The end will be one past the last
    auto current_execution_list_end = _perf_graph._execution_list_end.load(std::memory_order_relaxed);


    // The last entry in the current execution list
    auto current_execution_list_last = current_execution_list_end - 1 >= 0 ? current_execution_list_end - 1 : MAX_EXECUTION_LIST_SIZE;

//      std::cout << "current_execution_list_begin: " << current_execution_list_begin << std::endl;
//      std::cout << "current_execution_list_end: " << current_execution_list_end << std::endl;

    // This will synchronize with the thread_fence in addToExecutionList() so that all of the below
    // reads, will be reading synchronized memory
    std::atomic_thread_fence(std::memory_order_acquire);

    // Only happens if nothing has been added
    if (current_execution_list_end == 0 && _last_execution_list_end == current_execution_list_end)
      continue;

    // Iterate from the last thing printed (begin) to the last thing in the list (end)
    // If the time or memory of any section is above the threshold, print everything inbetween and
    // update begin

    // Are we still sitting in the same place as the last iteration?  If so, we need to print progress
    if (_last_execution_list_end == current_execution_list_end)
    {
      auto & last_section_increment = _execution_list[current_execution_list_last];

      // Is the last section to run still running?
      if (last_section_increment._state == PerfGraph::IncrementState::started)
      {
        if (!_printed_name_of_current_section)
        {
          // We need to print out everything on the stack before this that hasn't already been printed...
          for (unsigned int s = _printed_stack_level_end; s < _stack_level - 1; s++)
          {
            _console << std::string(s * 2, ' ') << _id_to_section_info[_print_thread_stack[s]._id]._live_message << '\n';
            _printed_stack_level_end++;
          }

          _console << std::string(2 * (_stack_level -1), ' ') << _id_to_section_info[last_section_increment._id]._live_message;

          _printed_name_of_current_section = true;

          _printed_stack_level_end++;
        }
        else // Need to print dots
          if (_id_to_section_info[last_section_increment._id]._print_dots)
            _console << " .";
      }
      else // If it's not, then we need to continue the section _before_ this
      {
        if (_stack_level > 0)
        {
          auto & last_stack_section = _print_thread_stack[_stack_level - 1];

          if (!_printed_name_of_current_section)
          {
            _console << "Still " << _id_to_section_info[last_stack_section._id]._live_message;

            _printed_name_of_current_section = true;
          }
          else // Need to print dots
            if (_id_to_section_info[last_stack_section._id]._print_dots)
              _console << " .";
        }
      }
    }

    // Current position in the execution list
    auto p = _last_execution_list_end;

//      std::cout << "p: " << p << std::endl;

    while (p != current_execution_list_end /*(current_execution_list_end + 1 < MAX_EXECUTION_LIST_SIZE ? current_execution_list_end + 1 : 0)*/)
    {
//        std::cout << "p: " << p << std::endl;

      auto next_p = p + 1 < MAX_EXECUTION_LIST_SIZE ? p + 1 : 0;
//        std::cout << "p: " << p << std::endl;
      auto & section_increment = _execution_list[p];

//        std::cout << "p: " << p << " Looking at: " << id_to_section_name[section_increment._id] << " " << section_increment._state << std::endl;

      // Do we need to increase the stack?
      if (section_increment._state == PerfGraph::IncrementState::started)
      {
        // Store this increment in the stack
        _print_thread_stack[_stack_level] = section_increment;

        _stack_level++;
      }
      else // See if we need to print it
      {
        // Get the beginning information for this section... it is the thing currently on the top of the stack
        auto & section_increment_start = _print_thread_stack[_stack_level - 1];

        auto time_increment = std::chrono::duration<double>(section_increment._time - section_increment_start._time).count();
        auto memory_increment = section_increment._memory - section_increment_start._memory;

        // Do they trigger the criteria or is it something that we were already partially printing?
        if (time_increment > 1. || memory_increment > 100/* || printed_name_of_current_section*/)
        {
          // We need to print out everything on the stack before this that hasn't already been printed...
          for (unsigned int s = _printed_stack_level_end; s < _stack_level - 1; s++)
          {
            // Before printing more things - need to stop printing what we were printing before
            if (_printed_name_of_current_section)
            {
              _console << '\n';
              _printed_name_of_current_section = false;
            }

            _console << std::string(s * 2, ' ') << _id_to_section_info[_print_thread_stack[s]._id]._live_message << '\n';
            _printed_stack_level_end++;
          }

          // Now print this thing
          if (!_printed_name_of_current_section || !_id_to_section_info[section_increment._id]._print_dots)
            _console << std::string(2 * (_stack_level - 1), ' ') << "Finishing: " << _id_to_section_info[section_increment._id]._live_message;

          _console << ": " << time_increment << ", " << memory_increment << '\n';

          _printed_name_of_current_section = false;
        }

        _stack_level--;

        _printed_stack_level_end = std::min(_printed_stack_level_end, _stack_level);
      }

      p = next_p;
    }

    // Make sure that everything comes out on the console
    _console << std::flush;

    _last_execution_list_end = current_execution_list_end;
  }
}
