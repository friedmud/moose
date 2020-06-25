//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PerfGraph.h"

// MOOSE Includes
#include "PerfGuard.h"
#include "MooseError.h"

// Note: do everything we can to make sure this only gets #included
// in the .C file... this is a heavily templated header that we
// don't want to expose to EVERY file in MOOSE...
#include "VariadicTable.h"

// libMesh Includes
#include "libmesh/auto_ptr.h"

// System Includes
#include <chrono>

const std::string PerfGraph::ROOT_NAME = "Root";

PerfGraph::PerfGraph(const std::string & root_name, MooseApp & app)
    : ConsoleStreamInterface(app), _root_name(root_name), _current_position(0),_execution_list_begin(0), _execution_list_end(0),  _active(true)
{
  // Start the printing thread
  _print_thread = std::thread([this] {
    auto & console = this->_console;

    auto & execution_list = this->_execution_list;

    auto & id_to_section_name = this->_id_to_section_name;

    // The current nesting level
    unsigned int current_level = 0;

    unsigned int horizontal_position = 0;

    unsigned int last_execution_list_end = 0;

    unsigned int last_printed_level = 0;

    // This is going to hold the entries to be printed going backward from the
    // current entry to be printed
    std::array<unsigned int, MAX_STACK_SIZE> back_buffer;

    unsigned int back_buffer_position = 0;

    while(true)
    {
      std::this_thread::sleep_for(std::chrono::seconds(1));

      auto current_execution_list_begin = this->_execution_list_begin.load(std::memory_order_relaxed);

      auto current_execution_list_end = this->_execution_list_end.load(std::memory_order_relaxed);

      // The last entry in the current execution list
      auto current_execution_list_last = current_execution_list_end - 1 >= 0 ? current_execution_list_end - 1 : MAX_EXECUTION_LIST_SIZE;

//      std::cout << "current_execution_list_begin: " << current_execution_list_begin << std::endl;
//      std::cout << "current_execution_list_end: " << current_execution_list_end << std::endl;

      // This will synchronize with the thread_fence in addToExecutionList() so that all of the below
      // reads, will be reading synchronized memory
      std::atomic_thread_fence(std::memory_order_acquire);

      // Iterate from the last thing printed (begin) to the last thing in the list (end)
      // If the time or memory of any section is above the threshold, print everything inbetween and
      // update begin

      // Current position in the execution list
      int p = current_execution_list_begin;

      // Need to finish printing
      if (current_execution_list_begin != current_execution_list_end && horizontal_position > 0)
      {
        console << " " << execution_list[current_execution_list_begin]._time << ", " << execution_list[current_execution_list_begin]._memory << std::endl;
        horizontal_position = 0;
      }


      while (p != current_execution_list_end /*(current_execution_list_end + 1 < MAX_EXECUTION_LIST_SIZE ? current_execution_list_end + 1 : 0)*/)
      {
        auto next_p = p + 1 < MAX_EXECUTION_LIST_SIZE ? p + 1 : 0;
//        std::cout << "p: " << p << std::endl;
        auto & section_increment = execution_list[p];

        if ((next_p == current_execution_list_end && current_execution_list_end == last_execution_list_end) || (section_increment._state == IncrementState::finished && (section_increment._time > 1. || section_increment._memory > 100)))
        {
          // We need to find the path from the thing we're currently at, back to the last level printed

          auto q = p;

          back_buffer_position = 0;

          int back_level = 0;

          int lowest_back_level = 0;

          while (q != current_execution_list_begin)
          {
//            std::cout << "q: " << q << std::endl;

            auto next_pos = q - 1 >= 0 ? q - 1 : MAX_EXECUTION_LIST_SIZE;

            // This will only happen if the previous thing that ended was inside this scope (i.e. one level higher)
            if (execution_list[q]._state == IncrementState::finished && execution_list[next_pos]._state == IncrementState::finished)
              back_level++;
            else if (execution_list[q]._state == IncrementState::started && execution_list[next_pos]._state == IncrementState::started)
              back_level--;

            // If we actually moved up in levels, then we need to print out the thing at next_pos
            if (back_level < lowest_back_level)
            {
              lowest_back_level = back_level;

              back_buffer[back_buffer_position] = next_pos;

              back_buffer_position++;
            }

            q = next_pos;

//            std::cout << "next q: " << q << std::endl;
          }

          // Now print out everything in the back buffer
          for (int i = back_buffer_position - 1; i >= 0; i--)
          {
            auto & back_section_increment = execution_list[back_buffer[i]];

            console /*<< "pre-"*/ << id_to_section_name[back_section_increment._id] << ": " << back_section_increment._time << ", " << back_section_increment._memory << std::endl;
            horizontal_position = 0;
          }

          // Now print out the one that is over the limits
          if (section_increment._state == IncrementState::finished)
          {
            console << id_to_section_name[section_increment._id] << ": " << section_increment._time << ", " << section_increment._memory << std::endl;

            horizontal_position = 0;
          }

//          if (p != (current_execution_list_end + 1 < MAX_EXECUTION_LIST_SIZE ? current_execution_list_end + 1 : 0))
          current_execution_list_begin = next_p;
        }

        p = next_p;

//        std::cout << "new p: " << p << std::endl;
//        std::cout << "stop p: " << (current_execution_list_end + 1 < MAX_EXECUTION_LIST_SIZE ? current_execution_list_end + 1 : 0) << std::endl;
      }

      if (current_execution_list_end == last_execution_list_end)
      {
        auto & section_increment = execution_list[current_execution_list_last];

        // If we haven't printed anything yet, need to print the name
        if (horizontal_position == 0)
        {
          if (section_increment._state == IncrementState::started)
            console << id_to_section_name[section_increment._id] << std::flush;
          else // This will happen if the enclosing scope still has more work to do after an enclosed scope has finished
            console << "Still" << std::flush;

          horizontal_position++;
        }
        else // Already printed the name - so print a dot
        {
          console << " ." << std::flush;
          horizontal_position++;
        }
      }
      else
        horizontal_position = 0;

      last_execution_list_end = current_execution_list_end;

      this->_execution_list_begin.store(current_execution_list_begin, std::memory_order_relaxed);
    }
  });

  // Not done in the initialization list on purpose because this object needs to be complete first
  _root_node = libmesh_make_unique<PerfNode>(registerSection(ROOT_NAME, 0));

  MemoryUtils::Stats stats;
  MemoryUtils::getMemoryStats(stats);
  auto start_memory =
      MemoryUtils::convertBytes(stats._physical_memory, MemoryUtils::MemUnits::Megabytes);

  // Set the initial time
  _root_node->setStartTimeAndMemory(std::chrono::steady_clock::now(), start_memory);

  // Add a call
  _root_node->incrementNumCalls();

  _stack[0] = _root_node.get();
}

unsigned int
PerfGraph::registerSection(const std::string & section_name, unsigned int level)
{
  auto it = _section_name_to_id.lower_bound(section_name);

  // Is it already registered?
  if (it != _section_name_to_id.end() && it->first == section_name)
    return it->second;

  // It's not...
  auto id = _section_name_to_id.size();
  _section_name_to_id.emplace_hint(it, section_name, id);
  _id_to_section_name[id] = section_name;
  _id_to_level[id] = level;

  return id;
}

const std::string &
PerfGraph::sectionName(const PerfID id) const
{
  auto find_it = _id_to_section_name.find(id);

  if (find_it == _id_to_section_name.end())
    mooseError("PerfGraph cannot find a section name associated with id: ", id);

  return find_it->second;
}

unsigned long int
PerfGraph::getNumCalls(const std::string & section_name)
{
  updateTiming();

  auto section_it = _section_time.find(section_name);

  if (section_it == _section_time.end())
    mooseError(
        "Unknown section_name: ",
        section_name,
        " in PerfGraph::getNumCalls()\nIf you are attempting to retrieve the root use \"Root\".");

  return section_it->second._num_calls;
}

Real
PerfGraph::getTime(const TimeType type, const std::string & section_name)
{
  updateTiming();

  auto section_it = _section_time.find(section_name);

  if (section_it == _section_time.end())
    mooseError(
        "Unknown section_name: ",
        section_name,
        " in PerfGraph::getTime()\nIf you are attempting to retrieve the root use \"Root\".");

  auto app_time = _section_time_ptrs[0]->_total;

  switch (type)
  {
    case SELF:
      return section_it->second._self;
    case CHILDREN:
      return section_it->second._children;
    case TOTAL:
      return section_it->second._total;
    case SELF_AVG:
      return section_it->second._self / static_cast<Real>(section_it->second._num_calls);
    case CHILDREN_AVG:
      return section_it->second._children / static_cast<Real>(section_it->second._num_calls);
    case TOTAL_AVG:
      return section_it->second._total / static_cast<Real>(section_it->second._num_calls);
    case SELF_PERCENT:
      return 100. * (section_it->second._self / app_time);
    case CHILDREN_PERCENT:
      return 100. * (section_it->second._children / app_time);
    case TOTAL_PERCENT:
      return 100. * (section_it->second._total / app_time);
    case SELF_MEMORY:
      return section_it->second._self_memory;
    case CHILDREN_MEMORY:
      return section_it->second._children_memory;
    case TOTAL_MEMORY:
      return section_it->second._total_memory;
    default:
      ::mooseError("Unknown TimeType");
  }
}

void
PerfGraph::push(const PerfID id)
{
  if (!_active)
    return;

  auto new_node = _stack[_current_position]->getChild(id);

  MemoryUtils::Stats stats;
  MemoryUtils::getMemoryStats(stats);
  auto start_memory =
      MemoryUtils::convertBytes(stats._physical_memory, MemoryUtils::MemUnits::Megabytes);

  // Set the start time
  auto current_time = std::chrono::steady_clock::now();

  new_node->setStartTimeAndMemory(current_time, start_memory);

  // Increment the number of calls
  new_node->incrementNumCalls();

  _current_position++;

  if (_current_position >= MAX_STACK_SIZE)
    mooseError("PerfGraph is out of stack space!");

  _stack[_current_position] = new_node;

  // Add this to the exection list
  addToExecutionList(id, IncrementState::started, 0., 0.);
}

void
PerfGraph::pop()
{
  if (!_active)
    return;

  MemoryUtils::Stats stats;
  MemoryUtils::getMemoryStats(stats);
  auto now_memory =
      MemoryUtils::convertBytes(stats._physical_memory, MemoryUtils::MemUnits::Megabytes);

  auto current_time = std::chrono::steady_clock::now();

  auto & current_node = _stack[_current_position];

  auto time_increment =  std::chrono::duration<double>(current_time - current_node->startTime()).count();
  long int memory_increment = now_memory - current_node->startMemory();

  current_node->addTimeAndMemory(current_time, now_memory);

  _current_position--;

  // Add this to the exection list
  addToExecutionList(current_node->id(), IncrementState::finished, time_increment, memory_increment);
}

void
PerfGraph::updateTiming()
{
  // First update all of the currently running nodes
  auto now = std::chrono::steady_clock::now();

  MemoryUtils::Stats stats;
  MemoryUtils::getMemoryStats(stats);
  auto now_memory =
      MemoryUtils::convertBytes(stats._physical_memory, MemoryUtils::MemUnits::Megabytes);

  for (unsigned int i = 0; i <= _current_position; i++)
  {
    auto node = _stack[i];
    node->addTimeAndMemory(now, now_memory);
    node->setStartTimeAndMemory(now, now_memory);
  }

  // Zero out the entries
  for (auto & section_time_it : _section_time)
  {
    auto & section_time = section_time_it.second;

    section_time._num_calls = 0;
    section_time._self = 0.;
    section_time._children = 0.;
    section_time._total = 0.;
    section_time._self_memory = 0;
    section_time._children_memory = 0;
    section_time._total_memory = 0.;
  }

  recursivelyFillTime(_root_node.get());

  // Update vector pointing to section times
  // Note: we are doing this _after_ recursively filling
  // because new entries may have been created
  _section_time_ptrs.resize(_id_to_section_name.size());

  for (auto & section_time_it : _section_time)
  {
    auto id = _section_name_to_id[section_time_it.first];

    _section_time_ptrs[id] = &section_time_it.second;
  }
}

void
PerfGraph::recursivelyFillTime(PerfNode * current_node)
{
  auto id = current_node->id();

  auto self = std::chrono::duration<double>(current_node->selfTime()).count();
  auto children = std::chrono::duration<double>(current_node->childrenTime()).count();
  auto total = std::chrono::duration<double>(current_node->totalTime()).count();
  auto num_calls = current_node->numCalls();

  auto self_memory = current_node->selfMemory();
  auto children_memory = current_node->childrenMemory();
  auto total_memory = current_node->totalMemory();

  // RHS insertion on purpose
  auto & section_time = _section_time[_id_to_section_name[id]];

  section_time._self += self;
  section_time._children += children;
  section_time._total += total;
  section_time._num_calls += num_calls;

  section_time._self_memory += self_memory;
  section_time._children_memory += children_memory;
  section_time._total_memory += total_memory;

  for (auto & child_it : current_node->children())
    recursivelyFillTime(child_it.second.get());
}

void
PerfGraph::recursivelyPrintGraph(PerfNode * current_node,
                                 FullTable & vtable,
                                 unsigned int level,
                                 unsigned int current_depth)
{
  mooseAssert(_id_to_section_name.find(current_node->id()) != _id_to_section_name.end(),
              "Unable to find section name!");

  auto & name = current_node->id() == 0 ? _root_name : _id_to_section_name[current_node->id()];

  mooseAssert(_id_to_level.find(current_node->id()) != _id_to_level.end(), "Unable to find level!");
  auto & node_level = _id_to_level[current_node->id()];

  if (node_level <= level)
  {
    mooseAssert(!_section_time_ptrs.empty(),
                "updateTiming() must be run before recursivelyPrintGraph!");

    auto section = std::string(current_depth * 2, ' ') + name;

    // The total time of the root node
    auto total_root_time = _section_time_ptrs[0]->_total;

    auto num_calls = current_node->numCalls();
    auto self = std::chrono::duration<double>(current_node->selfTime()).count();
    auto self_avg = self / static_cast<Real>(num_calls);
    auto self_percent = 100. * self / total_root_time;

    auto children = std::chrono::duration<double>(current_node->childrenTime()).count();
    auto children_avg = children / static_cast<Real>(num_calls);
    auto children_percent = 100. * children / total_root_time;

    auto total = std::chrono::duration<double>(current_node->totalTime()).count();
    auto total_avg = total / static_cast<Real>(num_calls);
    auto total_percent = 100. * total / total_root_time;

    auto self_memory = current_node->selfMemory();
    auto total_memory = current_node->totalMemory();

    vtable.addRow(section,
                  num_calls,
                  self,
                  self_avg,
                  self_percent,
                  self_memory,
                  //                  children,
                  //                  children_avg,
                  //                  children_percent,
                  total,
                  total_avg,
                  total_percent,
                  total_memory);

    current_depth++;
  }

  for (auto & child_it : current_node->children())
    recursivelyPrintGraph(child_it.second.get(), vtable, level, current_depth);
}

void
PerfGraph::recursivelyPrintHeaviestGraph(PerfNode * current_node,
                                         FullTable & vtable,
                                         unsigned int current_depth)
{
  mooseAssert(!_section_time_ptrs.empty(),
              "updateTiming() must be run before recursivelyPrintGraph!");

  auto & name = current_node->id() == 0 ? _root_name : _id_to_section_name[current_node->id()];

  auto section = std::string(current_depth * 2, ' ') + name;

  // The total time of the root node
  auto total_root_time = _section_time_ptrs[0]->_total;

  auto num_calls = current_node->numCalls();
  auto self = std::chrono::duration<double>(current_node->selfTime()).count();
  auto self_avg = self / static_cast<Real>(num_calls);
  auto self_percent = 100. * self / total_root_time;

  auto children = std::chrono::duration<double>(current_node->childrenTime()).count();
  auto children_avg = children / static_cast<Real>(num_calls);
  auto children_percent = 100. * children / total_root_time;

  auto total = std::chrono::duration<double>(current_node->totalTime()).count();
  auto total_avg = total / static_cast<Real>(num_calls);
  auto total_percent = 100. * total / total_root_time;

  auto self_memory = current_node->selfMemory();
  auto total_memory = current_node->totalMemory();

  vtable.addRow(section,
                num_calls,
                self,
                self_avg,
                self_percent,
                self_memory,
                //                children,
                //                children_avg,
                //                children_percent,
                total,
                total_avg,
                total_percent,
                total_memory);

  current_depth++;

  if (!current_node->children().empty())
  {
    PerfNode * heaviest_child = nullptr;

    for (auto & child_it : current_node->children())
    {
      auto current_child = child_it.second.get();

      if (!heaviest_child || (current_child->totalTime() > heaviest_child->totalTime()))
        heaviest_child = current_child;
    }

    recursivelyPrintHeaviestGraph(heaviest_child, vtable, current_depth);
  }
}

void
PerfGraph::print(const ConsoleStream & console, unsigned int level)
{
  updateTiming();

  console << "\nPerformance Graph:\n";
  FullTable vtable({"Section",
                    "Calls",
                    "Self(s)",
                    "Avg(s)",
                    "%",
                    "Mem(MB)",
                    // "Children(s)",
                    // "Avg(s)",
                    // "%",
                    "Total(s)",
                    "Avg(s)",
                    "%",
                    "Mem(MB)"},
                   10);

  vtable.setColumnFormat({
      VariadicTableColumnFormat::AUTO,    // Section Name
      VariadicTableColumnFormat::AUTO,    // Calls
      VariadicTableColumnFormat::FIXED,   // Self
      VariadicTableColumnFormat::FIXED,   // Avg.
      VariadicTableColumnFormat::PERCENT, // %
      VariadicTableColumnFormat::AUTO,    // Memory
      // VariadicTableColumnFormat::FIXED,     // Children
      // VariadicTableColumnFormat::FIXED,     // Avg.
      // VariadicTableColumnFormat::PERCENT,   // %
      VariadicTableColumnFormat::FIXED,   // Total
      VariadicTableColumnFormat::FIXED,   // Avg.
      VariadicTableColumnFormat::PERCENT, // %
      VariadicTableColumnFormat::AUTO,    // Memory
  });                                     // %

  vtable.setColumnPrecision({
      1, // Section Name
      0, // Calls
      3, // Self
      3, // Avg.
      2, // %
      0, // Memory
      3, // Total
      3, // Avg.
      2, // %
      0, // Memory
  });

  recursivelyPrintGraph(_root_node.get(), vtable, level);
  vtable.print(console);
}

void
PerfGraph::printHeaviestBranch(const ConsoleStream & console)
{
  updateTiming();

  console << "\nHeaviest Branch:\n";
  FullTable vtable({"Section",
                    "Calls",
                    "Self(s)",
                    "Avg(s)",
                    "%",
                    "Mem(MB)",
                    // "Children(s)",
                    // "Avg(s)",
                    // "%",
                    "Total(s)",
                    "Avg(s)",
                    "%",
                    "Mem(MB)"},
                   10);

  vtable.setColumnFormat({VariadicTableColumnFormat::AUTO,    // Section Name
                          VariadicTableColumnFormat::AUTO,    // Calls
                          VariadicTableColumnFormat::FIXED,   // Self
                          VariadicTableColumnFormat::FIXED,   // Avg.
                          VariadicTableColumnFormat::PERCENT, // %
                          VariadicTableColumnFormat::AUTO,    // Memory
                          // VariadicTableColumnFormat::FIXED,     // Children
                          // VariadicTableColumnFormat::FIXED,     // Avg.
                          // VariadicTableColumnFormat::PERCENT,   // %
                          VariadicTableColumnFormat::FIXED,   // Total
                          VariadicTableColumnFormat::FIXED,   // Avg.
                          VariadicTableColumnFormat::PERCENT, // %
                          VariadicTableColumnFormat::AUTO});  // Memory

  vtable.setColumnPrecision({
      1, // Section Name
      0, // Calls
      3, // Self
      3, // Avg.
      2, // %
      0, // Memory
      3, // Total
      3, // Avg.
      2, // %
      0, // Memory
  });

  recursivelyPrintHeaviestGraph(_root_node.get(), vtable);
  vtable.print(console);
}

void
PerfGraph::printHeaviestSections(const ConsoleStream & console, const unsigned int num_sections)
{
  updateTiming();

  console << "\nHeaviest Sections:\n";

  // Indirect Sort The Self Time
  std::vector<size_t> sorted;
  Moose::indirectSort(_section_time_ptrs.begin(),
                      _section_time_ptrs.end(),
                      sorted,
                      [](SectionTime * lhs, SectionTime * rhs) {
                        if (lhs && rhs)
                          return lhs->_self > rhs->_self;

                        // If the LHS exists - it's definitely bigger than a non-existant RHS
                        if (lhs)
                          return true;

                        // Both don't exist - so it doesn't matter how we sort them
                        return false;
                      });

  HeaviestTable vtable({"Section", "Calls", "Self(s)", "Avg.", "%", "Mem(MB)"}, 10);

  vtable.setColumnFormat({VariadicTableColumnFormat::AUTO, // Doesn't matter
                          VariadicTableColumnFormat::AUTO,
                          VariadicTableColumnFormat::FIXED,
                          VariadicTableColumnFormat::FIXED,
                          VariadicTableColumnFormat::PERCENT,
                          VariadicTableColumnFormat::AUTO});

  vtable.setColumnPrecision({1, 1, 3, 3, 2, 1});

  mooseAssert(!_section_time_ptrs.empty(),
              "updateTiming() must be run before printHeaviestSections()!");

  // The total time of the root node
  auto total_root_time = _section_time_ptrs[0]->_total;

  // Now print out the largest ones
  for (unsigned int i = 0; i < num_sections; i++)
  {
    auto id = sorted[i];

    vtable.addRow(id == 0 ? _root_name : _id_to_section_name[id],
                  _section_time_ptrs[id]->_num_calls,
                  _section_time_ptrs[id]->_total_memory,
                  _section_time_ptrs[id]->_self,
                  _section_time_ptrs[id]->_self /
                      static_cast<Real>(_section_time_ptrs[id]->_num_calls),
                  100 * _section_time_ptrs[id]->_self / total_root_time);
  }

  vtable.print(console);
}
