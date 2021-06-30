//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

// Forward Declarations
class PerfGraph;
class PerfGraphLivePrint;
class PerfNode;

namespace moose
{
namespace internal
{

// Forward Declarations
class PerfGraphRegistry;

/**
 * Get the global PerfGraphRegistry singleton.
 */
PerfGraphRegistry & getPerfGraphRegistry();

/**
 * The place where all timed sections will be stored
 */
class PerfGraphRegistry
{
public:
  /**
   * Used to hold metadata about the registered sections
   */
  struct SectionInfo
  {
    /// Unique ID
    PerfID _id;

    /// The name
    std::string _name;

    /// Print level (verbosity level)
    unsigned int _level;

    /// Message to print while the section is running
    std::string _live_message;

    /// Whether or not to print dots while this section runs
    bool _print_dots;
  };

  /**
   * Call to register a named section for timing.
   *
   * @param section_name The name of the code section to be timed
   * @param level The importance of the timer - lower is more important (0 will always come out)
   * @return The ID of the section - use when starting timing
   */
  PerfID registerSection(const std::string & section_name, const unsigned int level);

  /**
   * Call to register a named section for timing.
   *
   * @param section_name The name of the code section to be timed
   * @param level The importance of the timer - lower is more important (0 will always come out)
   * @param live_message The message to be printed to the screen during execution
   * @param print_dots Whether or not progress dots should be printed for this section
   * @return The ID of the section - use when starting timing
   */
  PerfID registerSection(const std::string & section_name,
                         const unsigned int level,
                         const std::string & live_message,
                         const bool print_dots = true);

protected:
  PerfGraphRegistry(){};

  /**
   * The internal function that actually carries out the registration
   */
  PerfID actuallyRegisterSection(const std::string & section_name,
                                 const unsigned int level,
                                 const std::string & live_message,
                                 const bool print_dots = true);

  /// Map of section names to IDs
  std::unordered_map<std::string, PerfID> _section_name_to_id;

  /// Map of IDs to section information
  std::unordered_map<PerfID, SectionInfo> _id_to_section_info;

  /// So it can be constructed
  friend PerfGraphRegistry & getPerfGraphRegistry();

  /// For access from the PerfGraph system
  friend class ::PerfGraph;
  friend class ::PerfGraphLivePrint;
  friend class ::PerfNode;
};

}
}
