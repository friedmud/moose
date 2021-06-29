//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE Includes
#include "MooseTypes.h"
#include "PerfNode.h"
#include "IndirectSort.h"
#include "ConsoleStream.h"
#include "ConsoleStreamInterface.h"
#include "MooseError.h"
#include "MemoryUtils.h"
#include "PerfGraphRegistry.h"

// System Includes
#include <array>
#include <atomic>
#include <thread>
#include <future>
#include <mutex>

// Forward Declarations
class PerfGuard;
class PerfGraphLivePrint;

template <class... Ts>
class VariadicTable;

#define MAX_STACK_SIZE 100
#define MAX_EXECUTION_LIST_SIZE 10000

/**
 * The PerfGraph will hold the master list of all registered performance segments and
 * the head PerfNode
 */
class PerfGraph : protected ConsoleStreamInterface
{
public:
  using PerfGraphRegistry = moose::internal::PerfGraphRegistry;
  using SectionInfo = PerfGraphRegistry::SectionInfo;

  /**
   * For retrieving values
   */
  enum TimeType : char
  {
    SELF,
    CHILDREN,
    TOTAL,
    SELF_AVG,
    CHILDREN_AVG,
    TOTAL_AVG,
    SELF_PERCENT,
    CHILDREN_PERCENT,
    TOTAL_PERCENT,
    SELF_MEMORY,
    CHILDREN_MEMORY,
    TOTAL_MEMORY
  };

  /**
   * Create a new PerfGraph
   */
  PerfGraph(const std::string & root_name, MooseApp & app, const bool live_all);

  /**
   * Destructor
   */
  ~PerfGraph();

  /**
   * Print the tree out
   *
   * @param console The output stream to output to
   * @param level The log level, the higher the number the more output you get
   */
  void print(const ConsoleStream & console, unsigned int level);

  /**
   * Print out the heaviest branch through the tree
   *
   * @param console The output stream to output to
   */
  void printHeaviestBranch(const ConsoleStream & console);

  /**
   * Print out the heaviest sections that were timed
   *
   * @param console The output stream to output to
   */
  void printHeaviestSections(const ConsoleStream & console, const unsigned int num_sections);

  /**
   * Grab the name of a section
   */
  const std::string & sectionName(const PerfID id) const;

  /**
   * Whether or not timing is active
   *
   * When not active no timing information will be kept
   */
  bool active() const { return _active; }

  /**
   * Turn on or off timing
   */
  void setActive(bool active) { _active = active; }

  /**
   * Turn on or off live printing (if timing is off then live printing will be off too)
   */
  void setLivePrintActive(bool active) { _live_print_active = active; }

  /**
   * Forces all sections to be output live
   */
  void setLivePrintAll(bool active) { _live_print_all = active; }

  /**
   * Set the time limit before a message prints
   */
  void setLiveTimeLimit(Real time_limit){ _live_print_time_limit = time_limit; }

  /**
   * Sert the memory limit before a message prints
   */
  void setLiveMemoryLimit(unsigned int mem_limit) { _live_print_mem_limit = mem_limit; }

  /**
   * Get the number of calls for a section
   */
  unsigned long int getNumCalls(const std::string & section_name);

  /**
   * Get a reference to the time for a section
   */
  Real getTime(const TimeType type, const std::string & section_name);

  /**
   * Get a reference to the self time for a section
   *
   * This reference can be held onto and the value
   * will be updated anytime updateTiming() is called.
   */
  const Real & getSelfTime(const std::string & section_name)
  {
    return _section_time[section_name]._self;
  }

  /**
   * Get a reference to the children time for a section
   *
   * This reference can be held onto and the value
   * will be updated anytime updateTiming() is called.
   */
  const Real & getChildrenTime(const std::string & section_name)
  {
    return _section_time[section_name]._children;
  }

  /**
   * Get a reference to the total time for a section
   *
   * This reference can be held onto and the value
   * will be updated anytime updateTiming() is called.
   */
  const Real & getTotalTime(const std::string & section_name)
  {
    return _section_time[section_name]._total;
  }

  /**
   * Get a reference to the self memory usage for a section
   *
   * This reference can be held onto and the value
   * will be updated anyting updateTiming() is called
   */
  const long int & getSelfMemory(const std::string & section_name)
  {
    return _section_time[section_name]._self_memory;
  }

  /**
   * Get a reference to the children memory usage for a section
   *
   * This reference can be held onto and the value
   * will be updated anyting updateTiming() is called
   */
  const long int & getChildrenMemory(const std::string & section_name)
  {
    return _section_time[section_name]._children_memory;
  }

  /**
   * Get a reference to the total memory usage for a section
   *
   * This reference can be held onto and the value
   * will be updated anyting updateTiming() is called
   */
  const long int & getTotalMemory(const std::string & section_name)
  {
    return _section_time[section_name]._total_memory;
  }

  /**
   * Udates the time section_time and time for all currently running nodes
   */
  void updateTiming();

protected:
  typedef VariadicTable<std::string,
                        unsigned long int,
                        Real,
                        Real,
                        Real,
                        long int,
                        Real,
                        Real,
                        Real,
                        long int>
      FullTable;

  typedef VariadicTable<std::string, unsigned long int, Real, Real, Real, long int> HeaviestTable;

  /**
   * Use to hold the time for each section
   *
   * These will be filled by updateTiming()
   */
  struct SectionTime
  {
    Real _self = 0.;
    Real _children = 0.;
    Real _total = 0.;
    unsigned long int _num_calls = 0;
    long int _self_memory = 0;
    long int _children_memory = 0;
    long int _total_memory = 0;
  };

  // Whether or not an increment is the start of the finish increment
  enum IncrementState
  {
    started,
    printed,
    continued,
    finished
  };

  /**
   * Use to hold an increment of time and memory for a section
   * This is used in the TimedPrint capability.
   */
  struct SectionIncrement
  {
    PerfID _id;

    /// Whether or not this increment is the start of an increment or
    /// the finishing of an increment.
    IncrementState _state;
    unsigned int _print_stack_level;
    unsigned int _num_dots;
    std::chrono::time_point<std::chrono::steady_clock> _time;
    long int _memory;
  };

  /**
   * Add the information to the execution list
   *
   * Should only be called by push() and pop()
   */
  inline void addToExecutionList(const PerfID id,
                                 const IncrementState state,
                                 const std::chrono::time_point<std::chrono::steady_clock> time,
                                 const long int memory);

  /**
   * Add a Node onto the end of the end of the current callstack
   *
   * Note: only accessible by using PerfGuard!
   */
  void push(const PerfID id);

  /**
   * Remove a Node from the end of the current scope
   *
   * Note: only accessible by using PerfGuard!
   */
  void pop();

  /**
   * Helper for printing out the graph
   *
   * @param current_node The node to be working on right now
   * @param console Where to print to
   * @param level The level to print out below (<=)
   * @param current_depth - Used in the recursion
   */
  void recursivelyPrintGraph(PerfNode * current_node,
                             FullTable & vtable,
                             unsigned int level,
                             unsigned int current_depth = 0);

  /**
   * Helper for printing out the trace that has taken the most time
   *
   * @param current_node The node to be working on right now
   * @param console Where to print to
   * @param current_depth - Used in the recursion
   */
  void recursivelyPrintHeaviestGraph(PerfNode * current_node,
                                     FullTable & vtable,
                                     unsigned int current_depth = 0);

  /**
   * Updates the cumulative self/children/total time
   *
   * Note: requires that self/children/total time are resized and zeroed before calling.
   *
   * @param current_node The current node to work on
   */
  void recursivelyFillTime(PerfNode * current_node);

  /**
   * Helper for printing out the heaviest sections
   *
   * @param console Where to print to
   */
  void printHeaviestSections(const ConsoleStream & console);

  /// Whether or not to put everything in the perf graph
  bool _live_print_all;

  /// The PerfGraphRegistry
  PerfGraphRegistry & _perf_graph_registry;

  /// This processor id
  processor_id_type _pid;

  /// The name (handle) of the root node
  static const std::string ROOT_NAME;

  /// The name to display for the root of the graph when printing the entire graph
  const std::string _root_name;

  /// The root node of the graph
  std::unique_ptr<PerfNode> _root_node;

  /// The id for the root node
  PerfID _root_node_id;

  /// The current node position in the stack
  int _current_position;

  /// The full callstack.  Currently capped at a depth of 100
  std::array<PerfNode *, MAX_STACK_SIZE> _stack;

  /// A circular buffer for holding the execution list, this is read by the printing loop
  std::array<SectionIncrement, MAX_EXECUTION_LIST_SIZE> _execution_list;

  /// Where the print thread should start reading the execution list
  std::atomic<unsigned int> _execution_list_begin;

  /// Where the print thread should stop reading the execution list
  std::atomic<unsigned int> _execution_list_end;

  /// Map of section names to IDs
  std::unordered_map<std::string, PerfID> & _section_name_to_id;

  /// Map of IDs to section information
  std::unordered_map<PerfID, SectionInfo> & _id_to_section_info;

  /// The time for each section.  This is updated on updateTiming()
  /// Note that this is _total_ cumulative time across every place
  /// that section is in the graph
  ///
  /// I'm making this a map so that we can give out references to the values
  /// The three values are: self, children
  /// The map is on std::string because we might need to be able to retrieve
  /// timing values in a "late binding" situation _before_ the section
  /// has been registered.
  std::unordered_map<std::string, SectionTime> _section_time;

  /// Pointers into _section_time indexed on PerfID
  /// This is here for convenience and speed so we don't need
  /// to iterate over the above map much - and it makes it
  /// easier to sort
  std::vector<SectionTime *> _section_time_ptrs;

  /// Whether or not timing is active
  bool _active;

  /// Whether or not live printing is active
  bool _live_print_active;

  /// The promise to the print thread that will signal when to stop
  std::promise<bool> _done;

  /// Tell the print thread to teardown
  std::atomic<bool> _destructing;

  /// The object that is doing live printing
  std::unique_ptr<PerfGraphLivePrint> _live_print;

  /// The thread for printing sections as they execute
  std::thread _print_thread;

  /// The mutex to use with a condition_variable for waking up the print thread
  std::mutex _print_thread_mutex;

  /// The condition_variable to wake the print thread
  std::condition_variable _finished_section;

  /// The time limit before a message is printed (in seconds)
  Real _live_print_time_limit;

  /// The memory limit before a message is printed (in MB)
  unsigned int _live_print_mem_limit;

  // Here so PerfGuard is the only thing that can call push/pop
  friend class PerfGuard;
  friend class PerfGraphLivePrint;
};
