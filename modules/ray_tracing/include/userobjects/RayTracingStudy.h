#ifndef RAYTRACINGSTUDY_H
#define RAYTRACINGSTUDY_H

// Local Includes
#include "Ray.h"
#include "RayTracing.h"
#include "SendBuffer.h"
#include "ReceiveBuffer.h"
#include "RayTracingMethod.h"

// MOOSE Includes
#include "GeneralUserObject.h"
#include "LIFOBuffer.h"
#include "CircularBuffer.h"

// libMesh Includes
#include "libmesh/libmesh.h"
#include "libmesh/parallel.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h"

// System Includes
#include <memory>
#include <chrono>
#include <thread>
#include <stdio.h>
#include <algorithm>

// Forward Declarations
class FEProblem;
class RayKernel;
class RayProblemBase;
class RayTracingStudy;

template <>
InputParameters validParams<RayTracingStudy>();

/**
 * Base class for Ray tracing studies that will generate Rays and then propoagate
 * all of them to termination.
 *
 * Subclasses _must_ override generateRays()
 */
class RayTracingStudy : public GeneralUserObject
{
public:
  RayTracingStudy(const InputParameters & parameters);

  /**
   * A place to reset the starting point of rays before tracing
   */
  virtual void resetRays() {}

  /**
   * Generates Rays from each element in random directions and traces them out across the domain
   */
  void executeStudy();

  /**
   * Total number of processor crossings for Rays that finished on this processor
   */
  unsigned long int totalProcessorCrossings() { return _total_processor_crossings; }

  /**
   * Max number of processor crossings for Rays that finished on this processor
   */
  unsigned long int maxProcessorCrossings() { return _max_processor_crossings; }

  /**
   * Number of Ray/element intersections done by this processor
   */
  unsigned long int localIntersections() { return _local_intersections; }

  /**
   * Number of Ray/element intersections for rays that end on this processor
   */
  unsigned long int totalIntersections() { return _total_intersections; }

  /**
   * Total number of Ray/element integrations for rays that end on this processor
   */
  unsigned long int totalIntegrations() { return _total_integrations; }

  /**
   * Total amount of distance traveled by the rays
   */
  unsigned long int totalDistance() { return _total_distance; }

  /**
   * Total amount of distance integrated over by the rays
   */
  unsigned long int totalIntegratedDistance() { return _total_integrated_distance; }

  /**
   * Local Rays started
   */
  unsigned long int localRaysStarted() { return _rays_started; }

  /**
   * Total number of Rays started
   */
  unsigned long int totalRaysStarted() { return _all_rays_started; }

  /**
   * Rays traced on this processor
   */
  unsigned long int raysTraced() { return _rays_traced; }

  /**
   * Total number of received Ray buffers traced by this processor
   */
  unsigned long int rayBuffersReceived() { return _receive_buffer.buffersReceived(); }

  /**
   * Total number of times that messages were probed for
   */
  unsigned long int numProbes() { return _receive_buffer.numProbes(); }

  /**
   * Total number of rays that got fast laned
   */
  unsigned long int fastLaneRays() { return _receive_buffer.fastLaneRays(); }

  /**
   * Total number of ray buffers created in the ReceiveBuffer
   */
  unsigned long int receiveRayPoolCreated() { return _receive_buffer.rayPoolCreated(); }

  /**
   * Total number of buffers created in the ReceiveBuffer
   */
  unsigned long int receiveBufferPoolCreated() { return _receive_buffer.bufferPoolCreated(); }

  /**
   * Total number of buffers created in the SendBuffers
   */
  unsigned long int sendBufferPoolCreated()
  {
    unsigned long int total = 0;

    for (auto & buffer : _send_buffers)
      total += buffer.second->bufferPoolCreated();

    return total;
  }

  /**
   * Rays sent from this processor
   */
  unsigned long int raysSent()
  {
    unsigned long int total_sent = 0;

    for (auto & buffer : _send_buffers)
      total_sent += buffer.second->raysSent();

    return total_sent;
  }

  /**
   * Buffers sent from this processor
   */
  unsigned long int buffersSent()
  {
    unsigned long int total_sent = 0;

    for (auto & buffer : _send_buffers)
      total_sent += buffer.second->buffersSent();

    return total_sent;
  }

  /**
   * Total number of received Rays by this processor
   */
  unsigned long int raysReceived() { return _receive_buffer.raysReceived(); }

  /**
   * Total number of chunks traced
   */
  unsigned long int chunksTraced() { return _chunks_traced; }

  /**
   * Max distance of any ray
   */
  Real rayMaxDistance() { return _ray_max_distance; }

  /**
   * Set the maximum distance of any ray
   */
  void setRayMaxDistance(Real ray_max_distance) { _ray_max_distance = ray_max_distance; }

  /**
   * Get the total number of rays to run on
   */
  unsigned long int totalRays() { return _total_rays; }

  /**
   * Set the total number of rays
   */
  void setTotalRays(unsigned long int total_rays) { _total_rays = total_rays; }

  /**
   * Set the max ray distance
   */
  void setRayMaxDistance(unsigned long int ray_max_distance) { _ray_max_distance = ray_max_distance; }

  /**
   * Duration for execute() in seconds
   */
  Real executionTime() { return std::chrono::duration<Real>(_execution_time).count(); }

  /**
   * Duration for execute() in nanoseconds
   */
  Real executionTimeNano()
  {
    return std::chrono::duration<Real, std::nano>(_execution_time).count();
  }

  /**
   * Duration for creation of all Rays in seconds
   */
  Real generationTime() { return std::chrono::duration<Real>(_generation_time).count(); }

  /**
   * Duration for creation of all Rays in seconds
   */
  Real propagationTime() { return std::chrono::duration<Real>(_propagation_time).count(); }

  /**
   * Duration for tracing of all Rays in seconds
   */
  Real tracingTime() { return std::chrono::duration<Real>(_tracing_time).count(); }

  /**
   * Duration for buffering (and sending) of all Rays in seconds
   */
  Real bufferTime() { return std::chrono::duration<Real>(_buffer_time).count(); }

  /**
   * Duration for forced buffering (and sending) of all Rays in seconds
   */
  Real forceBufferTime() { return std::chrono::duration<Real>(_force_buffer_time).count(); }

  /**
   * Duration for receiving of all Rays in seconds
   */
  Real receiveTime() { return std::chrono::duration<Real>(_receive_time).count(); }

  /**
   * Amount of time (in seconds) spent in the loop checking for messages and creating Receives
   */
  Real receiveLoopTime() { return _receive_buffer.receiveLoopTime(); }

  /**
   * Amount of time (in seconds) spent finishing receives and collecting rays
   */
  Real cleanupRequestsTime() { return _receive_buffer.cleanupRequestsTime(); }

  /**
   * Time just spent busy waiting
   */
  Real waitTime() { return std::chrono::duration<Real>(_wait_time).count(); }

  /**
   * Time just spent communicating with the root processor
   */
  Real rootCommTime() { return std::chrono::duration<Real>(_root_comm_time).count(); }

  /**
   * Time spent by the root processor doing root processor things
   */
  Real rootTime() { return std::chrono::duration<Real>(_root_time).count(); }

  /**
   * Number of halo layers
   */
  unsigned long int haloSize() { return _halo_size; }

  /**
   * Number of Rays to buffer before flushing
   */
  unsigned long int maxBufferSize() { return _max_buffer_size; }

  /**
   * Set Number of Rays to buffer before flushing
   */
  void setMaxBufferSize(unsigned long int max_buffer_size) { _max_buffer_size = max_buffer_size; }

  /**
   * Number of Rays to generate before tracing
   */
  unsigned long int chunkSize() { return _chunk_size; }

  /**
   * Number of Rays to generate before tracing
   */
  void setChunkSize(unsigned long int chunk_size) { _chunk_size = chunk_size; }

  /**
   * Number of iterations before focing communication.
   */
  unsigned long int clicks() { return _clicks_per_communication; }

  /**
   * Set Number of iterations before forcing communication.
   */
  void setClicks(unsigned long int clicks_per_communication)
  {
    _clicks_per_communication = clicks_per_communication;
  }

  /**
   * Number of iterations before forcing root communication.
   */
  unsigned long int rootClicks() { return _clicks_per_root_communication; }

  /**
   * Set Number of iterations before forcing root communication.
   */
  void setRootClicks(unsigned long int clicks_per_root_communication)
  {
    _clicks_per_root_communication = clicks_per_root_communication;
  }

  /**
   * Number of iterations before looking for new rays
   */
  unsigned long int receiveClicks() { return _clicks_per_receive; }

  /**
   * Set Number of iterations before looking for new rays
   */
  void setReceiveClicks(unsigned long int clicks_per_receive)
  {
    _clicks_per_receive = clicks_per_receive;
    _receive_buffer.setReceiveClicks(clicks_per_receive);
  }

  /**
   * Set the ray tracing algorithm to use
   */
  void setMethod(RayTracingMethod method);

  /**
   * Get the current ray tracing algorithm
   */
  RayTracingMethod getMethod() { return _method; }

protected:
  /**
   * Subclasses should override this to determine how to generate Rays.
   * This will be called within execute() and makes up the "generation phase"
   * of the algorithm.
   */
  virtual void generateRays() = 0;

  /**
   * Second phase of Ray tracing.  After generating Rays this function
   * moves all the existing Rays through the domain.
   */
  void propagateRays();

  /**
   * Trace each Ray in the vector of Rays across the domain.
   * The Ray will either be buffered to send to the next processor
   * or it will be terminated if it reaches a boundary
   */
  void traceAndBuffer(std::vector<std::shared_ptr<Ray>>::iterator begin,
                      std::vector<std::shared_ptr<Ray>>::iterator end);

  /**
   * Traces out the given Ray's in chunks
   */
  void chunkyTraceAndBuffer(bool start_receives_only = false);

  /**
   * Reeive packets of Rays from other processors and trace
   * them out across the domain.
   */
  bool receiveAndTrace();

  /**
   * Flushes all Rays out of the send buffers
   */
  void flushBuffers();

  /**
   * Look for incoming messages about the number of rays finished
   */
  void receiveRaysFinished();

  /// The RayProblem
  RayProblemBase & _ray_problem;

  /// Number of groups
  unsigned long int _num_groups;

  /// The Mesh
  MooseMesh & _mesh;

  /// The Communicator
  const Parallel::Communicator & _comm;

  /// Size of the "halo" zone
  unsigned long int _halo_size;

  /// Number of rays to generate before tracing
  unsigned long int _chunk_size;

  /// Number of Rays to start in the volume
  unsigned long int _total_rays;

  /// Total volume of the domain
  Real _volume;

  /// Whether or not our working buffer is a LIFOBuffer
  bool _is_lifo;

  /// The set of rays we're currently working on
  WorkBufferType _working_buffer;

  /// Rays that are headed toward the center of the domain.
  /// They should be processed first
  WorkBufferType _fast_lane;

  /// Whether or not to use the fast lane buffer
  bool _use_fast_lane;

  /// Number of Rays to buffer before communication
  unsigned long int _max_buffer_size;

  /// Multiplier for the buffer size for growing the buffer
  Real _buffer_growth_multiplier;

  /// Multiplier for the buffer size for shrinking the buffer
  Real _buffer_shrink_multiplier;

  /// Minimum size of a SendBuffer
  unsigned long int _min_buffer_size;

  /// The rank of this processor (this actually takes time to lookup - so just do it once)
  processor_id_type _my_pid;

  /// Reusable Points
  Point _start;
  Point _end;
  Point _direction;

  /// Length of Rays
  Real _ray_length;

  /// Processor 0 will tell us when to stop
  unsigned long int _stop = 0;

  /// Number of Rays started on this processor
  unsigned long int _rays_started = 0;

  /// Number of Rays started everywhere
  unsigned long int _all_rays_started = 0;

  /// Number of Rays finished on this processor
  unsigned long int _rays_finished = 0;

  /// Temporary buffer used for transmission of _rays_finished
  unsigned long int _rays_finished_temp = 0;

  /// Number of Rays finished everywhere
  unsigned long int _all_rays_finished = 0;

  /// Total number of processor crossings for Rays that finished on this processor
  unsigned long int _total_processor_crossings = 0;

  /// Local number of Ray/element intersections
  unsigned long int _local_intersections = 0;

  /// Total number of Ray/element intersections
  unsigned long int _total_intersections = 0;

  /// Total number of Ray/element integrations
  unsigned long int _total_integrations = 0;

  /// Total distance traveled by all Rays
  Real _total_distance = 0;

  /// Total distance integrated over by all Rays
  Real _total_integrated_distance = 0;

  /// Max number of processor crossings for Rays that finished on this processor
  unsigned long int _max_processor_crossings = 0;

  /// Number of times a received Ray was traced
  unsigned long int _rays_traced = 0;

  /// Number of received Ray chunks traced over
  unsigned long int _chunks_traced = 0;

  /// How many processors have finished generating all of their rays
  processor_id_type _ranks_finished_generating = 0;

  std::unordered_map<processor_id_type, std::shared_ptr<SendBuffer>> _send_buffers;
  ReceiveBuffer _receive_buffer;

  /// Tags used in parallel communication
  Parallel::MessageTag _ray_buffer_tag = Parallel::MessageTag(RAY_BUFFER);

  /// Timing
  std::chrono::steady_clock::duration _execution_time = std::chrono::steady_clock::duration::zero();
  std::chrono::steady_clock::duration _generation_time =
      std::chrono::steady_clock::duration::zero();
  std::chrono::steady_clock::duration _propagation_time =
      std::chrono::steady_clock::duration::zero();
  std::chrono::steady_clock::duration _tracing_time = std::chrono::steady_clock::duration::zero();
  std::chrono::steady_clock::duration _buffer_time = std::chrono::steady_clock::duration::zero();
  std::chrono::steady_clock::duration _force_buffer_time =
      std::chrono::steady_clock::duration::zero();
  std::chrono::steady_clock::duration _receive_time = std::chrono::steady_clock::duration::zero();
  std::chrono::steady_clock::duration _wait_time = std::chrono::steady_clock::duration::zero();
  std::chrono::steady_clock::duration _root_comm_time = std::chrono::steady_clock::duration::zero();
  std::chrono::steady_clock::duration _root_time = std::chrono::steady_clock::duration::zero();

  /// Timouts for communication
  unsigned long int _clicks_per_communication;
  unsigned long int _clicks_per_root_communication;
  unsigned long int _clicks_per_receive;

  /// Max distance a Ray can travel before being killed
  Real _ray_max_distance;

  /// Whether or not to tolerate a Ray Tracing failure
  bool _tolerate_failure;

  /// Average angular flux for rays that have finished
  std::vector<Real> _average_finishing_angular_flux;

  /// Old Average angular flux for rays that have finished
  std::vector<Real> _old_average_finishing_angular_flux;

  /// The current ray tracing algorithm to use
  RayTracingMethod _method;

  ///// HARM Stuff ////

  /// Tag for sending rays finished
  Parallel::MessageTag _rays_finished_requests_tag = Parallel::MessageTag(RAYS_FINISHED);

  /// Whether or not rays_finished has been sent to each proc
  std::vector<bool> _rays_finished_requests_sent;

  /// Requests for sending the number of finished rays to every other processor
  std::vector<Parallel::Request> _rays_finished_requests;

  /// Values of rays destroyed on this processor that are being sent to other processors
  std::vector<unsigned long int> _rays_finished_requests_temps;

  /// Rays finished by each processor
  std::vector<unsigned long int> _rays_finished_per_proc;

private:
  /**
   * Propagate rays using SMART
   */
  void smartPropagate();

  /**
   * Propagate rays using HARM
   */
  void harmPropagate();

  /**
   * Propagate rays using BS
   */
  void bsPropagate();

  /**
   * These don't mean anythign here:
   */
  virtual void execute() {}
  virtual void initialize() {}
  virtual void finalize() {}

  /**
   * Timers
   */
  PerfID _execute_study_timer;
  PerfID _generate_timer;
  PerfID _propagate_timer;
  PerfID _trace_and_buffer_timer;
  PerfID _trace_timer;
  PerfID _buffer_timer;
};

#endif
