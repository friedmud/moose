#ifndef RANKMAP_H
#define RANKMAP_H

#include "PerfGraphInterface.h"

#include "libmesh/parallel_object.h"

/**
 * Builds lists and maps that help in knowing which physical hardware nodes each rank is on.
 *
 * Note: large chunks of this code were originally committed by @dschwen in PR #12351
 *
 * https://github.com/idaholab/moose/pull/12351
 */
class RankMap : ParallelObject, PerfGraphInterface
{
public:
  /**
   * Constructs and fills the map
   */
  RankMap(const Parallel::Communicator & comm, PerfGraph & perf_graph);

  /**
   * Returns the "hardware ID" (a unique ID given to each physical compute node in the job)
   * for a given processor ID (rank)
   */
  unsigned int hardwareID(processor_id_type pid) const { return _rank_to_hardware_id[pid]; }

  /**
   * Returns the ranks that are on the given hardwareID (phsical node in the job)
   */
  const std::vector<processor_id_type> & ranks(unsigned int hardware_id) const
  {
    return _hardware_id_to_ranks.at(hardware_id);
  }

  /**
   * Vector containing the hardware ID for each PID
   */
  const std::vector<unsigned int> & rankHardwareIds() const { return _rank_to_hardware_id; }

  /**
   * The local rank on the node
   */
  processor_id_type localRank() const { return _local_rank; }

protected:
  PerfID _construct_timer;

  /// Map of hardware_id -> ranks on that node
  std::map<unsigned int, std::vector<processor_id_type>> _hardware_id_to_ranks;

  /// Each entry corresponds to the hardware_id for that PID
  std::vector<unsigned int> _rank_to_hardware_id;

  /// The rank of the MPI process within its node
  processor_id_type _local_rank;
};

#endif
