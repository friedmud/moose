// Local Includes
#include "RayTracingStudy.h"
#include "RayProblem.h"
#include "TraceRay.h"

// MOOSE Includes
#include "ParallelUniqueId.h"
#include "MooseMesh.h"

// libMesh Includes
#include "libmesh/threads.h"

// System Includes
#include "unistd.h"
#include <numeric>

#define BATCH_SIZE 10

template <>
InputParameters
validParams<RayTracingStudy>()
{
  InputParameters params = validParams<GeneralUserObject>();

  params.addParam<unsigned int>(
      "halo_size", 1, "The size of the ghosted layer around each local set of elements");
  params.addRequiredParam<unsigned int>("num_rays", "The number of rays to use");
  params.addRequiredParam<Real>("ray_distance", "Rays will travel until they reach this distance");

  params.addParam<unsigned int>("send_buffer_size", 100, "The size of the send buffer");
  params.addParam<unsigned int>(
      "chunk_size", 100, "The number of rays to process at one time during generation");
  params.addParam<unsigned int>(
      "clicks_per_communication", 100, "Iterations to wait before communicating");
  params.addParam<unsigned int>(
      "clicks_per_root_communication", 10000, "Iterations to wait before communicating with root");
  params.addParam<unsigned int>(
      "clicks_per_receive", 1, "Iterations to wait before checking for new rays");

  params.addParam<bool>("blocking", false, "Whether to do blocking or non-blocking receives");

  params.addParam<bool>(
      "tolerate_failure", false, "Whether or not to tolerate a ray tracing failure");

  MooseEnum methods("smart harm bs", "smart");

  params.addParam<MooseEnum>("method", methods, "The ray tracing algorithm to use");

  return params;
}

RayTracingStudy::RayTracingStudy(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _ray_problem(dynamic_cast<RayProblemBase &>(*parameters.get<FEProblem *>("_fe_problem"))),
    _num_groups(_ray_problem.numGroups()),
    _mesh(_fe_problem.mesh()),
    _comm(_mesh.comm()),
    _halo_size(getParam<unsigned int>("halo_size")),
    _total_rays(getParam<unsigned int>("num_rays")),
    _max_buffer_size(getParam<unsigned int>("send_buffer_size")),
    _chunk_size(getParam<unsigned int>("chunk_size")),
    _my_pid(_comm.rank()),
    _ray_length(_ray_problem.domainMaxLength()),
    _receive_buffer(
        _comm, _mesh, getParam<unsigned int>("clicks_per_receive"), getParam<bool>("blocking")),
    _clicks_per_communication(getParam<unsigned int>("clicks_per_communication")),
    _clicks_per_root_communication(getParam<unsigned int>("clicks_per_root_communication")),
    _clicks_per_receive(getParam<unsigned int>("clicks_per_receive")),
    _ray_max_distance(getParam<Real>("ray_distance")),
    _tolerate_failure(getParam<bool>("tolerate_failure")),
    _average_finishing_angular_flux(_ray_problem.numGroups() * _ray_problem.numPolar()),
    _old_average_finishing_angular_flux(_ray_problem.numGroups() * _ray_problem.numPolar()),
    _method((RayTracingMethod)(int)(getParam<MooseEnum>("method")))
{
  setMethod(_method);
}

void
RayTracingStudy::executeStudy()
{
  // Reset the internal data

  /// Processor 0 will tell us when to stop
  _stop = 0;

  _rays_started = 0;
  _all_rays_started = 0;
  _rays_finished = 0;
  _all_rays_finished = 0;

  /// Total number of processor crossings for Rays that finished on this processor
  _total_processor_crossings = 0;

  /// Total number of Ray/element intersections
  _total_intersections = 0;

  /// Total number of Ray/element intersections
  _total_integrations = 0;

  /// Total distance traveled by all Rays
  _total_distance = 0;

  /// Total distance integrated over by all Rays
  _total_integrated_distance = 0;

  /// Max number of processor crossings for Rays that finished on this processor
  _max_processor_crossings = 0;

  /// Number of times a received Ray was traced
  _rays_traced = 0;

  /// Number of received Ray buffers traced over
  _ray_buffers_traced = 0;

  /// How many processors have finished generating all of their rays
  _ranks_finished_generating = 0;

  for (auto & buffer : _send_buffers)
    buffer.second->clear();

  _send_buffers.clear();
  _receive_buffer.clear();

  /// Timing
  _execution_time = std::chrono::steady_clock::duration::zero();
  _generation_time = std::chrono::steady_clock::duration::zero();
  _tracing_time = std::chrono::steady_clock::duration::zero();
  _buffer_time = std::chrono::steady_clock::duration::zero();
  _force_buffer_time = std::chrono::steady_clock::duration::zero();
  _receive_time = std::chrono::steady_clock::duration::zero();
  _wait_time = std::chrono::steady_clock::duration::zero();
  _root_comm_time = std::chrono::steady_clock::duration::zero();
  _root_time = std::chrono::steady_clock::duration::zero();

  /// Average angular flux for rays that have finished
  std::fill(_average_finishing_angular_flux.begin(), _average_finishing_angular_flux.end(), 0);

  /// Reset HARM stuff
  if (_method == HARM)
  {
    _rays_finished_requests.resize(_comm.size());

    _rays_finished_requests_sent.resize(_comm.size());
    std::fill(_rays_finished_requests_sent.begin(), _rays_finished_requests_sent.end(), false);

    _rays_finished_requests_temps.resize(_comm.size());
    std::fill(_rays_finished_requests_temps.begin(), _rays_finished_requests_temps.end(), 0);

    _rays_finished_per_proc.resize(_comm.size());
    std::fill(_rays_finished_per_proc.begin(), _rays_finished_per_proc.end(), 0);
  }

  /// Start Execution
  auto execution_start_time = std::chrono::steady_clock::now();
  ;

  // global_packing_time = std::chrono::steady_clock::duration::zero();
  // global_unpacking_time = std::chrono::steady_clock::duration::zero();

  //  _communicator.barrier();

  generateRays();

  //  _communicator.barrier();

  propagateRays();

  auto execution_end_time = std::chrono::steady_clock::now();

  _execution_time += execution_end_time - execution_start_time;

  /// Compute and keep the average finishing angular flux (to use as the starting values next time)
  _old_average_finishing_angular_flux = _average_finishing_angular_flux;
  _comm.sum(_old_average_finishing_angular_flux);
  for (auto & v : _old_average_finishing_angular_flux)
    v /= _all_rays_finished;
}

void
RayTracingStudy::setMethod(RayTracingMethod method)
{
  _method = method;

  _receive_buffer.setMethod(_method);

  for (auto & send_buffer_iter : _send_buffers)
    send_buffer_iter.second->setMethod(_method);
}

template <typename T>
void
nonblockingSum(const libMesh::Parallel::Communicator & comm,
               T & r,
               T & o,
               libMesh::Parallel::Request & req)
{
  if (comm.size() > 1)
    libmesh_call_mpi(MPI_Iallreduce(
        &r, &o, 1, libMesh::Parallel::StandardType<T>(&r), MPI_SUM, comm.get(), req.get()));
  else
    o = r;
}

void
RayTracingStudy::propagateRays()
{
  std::cout << "Method: " << _method << std::endl;

  switch (_method)
  {
    case SMART:
      smartPropagate();
      break;
    case HARM:
      harmPropagate();
      break;
    case BS:
      bsPropagate();
      break;
    default:
      mooseError("Unknown RayTracingMethod!");
  }
}

void
RayTracingStudy::traceAndBuffer(std::vector<std::shared_ptr<Ray>> & rays)
{
  auto num_rays = rays.size();

#pragma omp parallel for schedule(guided)
  for (auto ray_index = decltype(num_rays)(0); ray_index < num_rays; ray_index++)
  {
    auto & ray = rays[ray_index];

    if (ray->startingElem()->processor_id() == _my_pid) // This will be false if we've picked up a
                                                        // banked ray that needs to start on another
                                                        // processor
    {
      RayProblemTraceRay tr(_ray_problem,
                            _mesh.getMesh(),
                            _halo_size,
                            _ray_max_distance,
                            _ray_length,
                            _tolerate_failure,
                            omp_get_thread_num());

      tr.trace(ray);
    }
  }

  // Figure out where they went
  for (auto & ray : rays)
  {
    // Means that the Ray needs to be given to another processor to finish tracing
    if (ray->shouldContinue())
    {
      ray->processorCrossings()++;

      processor_id_type next_pid = ray->startingElem()->processor_id();

      if (!_send_buffers[next_pid])
        _send_buffers[next_pid] =
            std::make_shared<SendBuffer>(_comm, next_pid, _max_buffer_size, _method);

      _send_buffers[next_pid]->addRay(ray);
    }
    else
    {
      unsigned long int pc = ray->processorCrossings();

      _total_processor_crossings += pc;

      _max_processor_crossings = std::max(_max_processor_crossings, pc);

      _total_intersections += ray->intersections();

      _total_integrations +=
          (ray->intersections() * ray->polarSins().size() * _ray_problem.numGroups());

      _total_distance += ray->distance();

      _total_integrated_distance += ray->integratedDistance();

      for (unsigned int i = 0; i < _ray_problem.numGroups(); i++)
        _average_finishing_angular_flux[i] += ray->_data[i];

      _rays_finished++;
    }
  }

  if (_method == HARM)
    flushBuffers();

  // Try to finish off some sends
  for (auto & send_buffer_iter : _send_buffers)
    send_buffer_iter.second->cleanupRequests();
}

void
RayTracingStudy::chunkyTraceAndBuffer(std::vector<std::shared_ptr<Ray>> & rays)
{
  auto num_rays = rays.size();

  unsigned int num_chunks = num_rays / _chunk_size;

  if (_method == HARM || _method == BS || num_chunks == 0)
    num_chunks = 1;

  auto current_beginning = rays.begin();

  for (auto chunk = decltype(num_chunks)(0); chunk < num_chunks; chunk++)
  {
    auto current_chunk_size = _chunk_size;

    if (_chunk_size < num_rays)
    {
      // The last chunk gets the rest of the rays
      if (chunk == num_chunks - 1)
        current_chunk_size += num_rays % _chunk_size;
    }
    else
      current_chunk_size = num_rays;

    std::vector<std::shared_ptr<Ray>> current_rays(current_beginning,
                                                   current_beginning + current_chunk_size);

    current_beginning += current_chunk_size;

    traceAndBuffer(current_rays);

    // Interleave receiving and tracing of incoming rays
    // by doing this while we generate Rays we are overlapping communication and computation
    // This call makes this recursive...
    if (_method == SMART)
      receiveAndTrace();
  }
}

bool
RayTracingStudy::receiveAndTrace()
{
  bool traced_some = false;

  _receive_buffer.receive();

  unsigned int received_count = 0;

  while (true)
  {
    auto rays = _receive_buffer.getRays();

    // Make sure we have some Rays to trace!
    if (!rays)
      break;

    traced_some = true;

    received_count += rays->size();

    _rays_traced += rays->size();
    _ray_buffers_traced++;

    if (_method == SMART)
      chunkyTraceAndBuffer(*rays);
    else
      traceAndBuffer(*rays);
  }

  return traced_some;
}

void
RayTracingStudy::flushBuffers()
{
  for (auto & send_buffer_iter : _send_buffers)
    send_buffer_iter.second->forceSend();
}

void
RayTracingStudy::smartPropagate()
{
  // std::cout << "Using SMART!!!!" << std::endl;

  Parallel::Request rays_started_request;
  Parallel::Request rays_finished_request;

  // Get the number of rays that were started in the whole domain
  nonblockingSum(_comm, _rays_started, _all_rays_started, rays_started_request);

  // Good time to get rid of whatever's currently in our SendBuffers
  flushBuffers();

  // Use these to try to delay some forced communication
  unsigned int non_tracing_clicks = 0;

  bool traced_some = true;

  // Keep bouncing the rays around until they've all completed
  while (true)
  {
    traced_some = receiveAndTrace();

    bool sending = false;

    for (auto & send_buffer : _send_buffers)
      sending = sending || send_buffer.second->currentlySending() ||
                send_buffer.second->currentlyBuffered();

    if (traced_some || _receive_buffer.currentlyReceiving())
      non_tracing_clicks = 0;
    else
      non_tracing_clicks++;

    if (_method == SMART && non_tracing_clicks >= _clicks_per_communication)
    {
      non_tracing_clicks = 0;
      flushBuffers();
    }

    if (!(traced_some || _receive_buffer.currentlyReceiving() || sending) &&
        rays_started_request.test() && rays_finished_request.test())
    {
      if (_all_rays_started == _all_rays_finished)
        return;

      _rays_finished_temp = _rays_finished;
      nonblockingSum(_comm, _rays_finished_temp, _all_rays_finished, rays_finished_request);
    }
  }
}

void
RayTracingStudy::harmPropagate()
{
  // std::cout << "Using Harm!!!" << std::endl;

  Parallel::Request rays_started_request;
  Parallel::Request rays_finished_request;

  // Get the number of rays that were started in the whole domain
  nonblockingSum(_comm, _rays_started, _all_rays_started, rays_started_request);

  // All rays have been traced out, so time to communicate
  flushBuffers();

  // HARM only does some communication based on times through the loop.
  // This counter will be used for that
  unsigned int communication_clicks;

  // Cache these
  processor_id_type num_procs = _comm.size();
  processor_id_type pid = _comm.rank();

  Parallel::Status rays_finished_probe_status;
  int rays_finished_probe_flag;

  // Keep bouncing the rays around until they've all completed
  while (true)
  {
    receiveAndTrace();

    flushBuffers();

    if (communication_clicks > num_procs)
    {
      // Receive messages about rays being finished
      do
      {
        MPI_Iprobe(MPI_ANY_SOURCE,
                   _rays_finished_requests_tag.value(),
                   _comm.get(),
                   &rays_finished_probe_flag,
                   rays_finished_probe_status.get());

        if (rays_finished_probe_flag)
        {
          auto proc = rays_finished_probe_status.source();
          _comm.receive(proc, _rays_finished_per_proc[proc], _rays_finished_requests_tag);
        }
      } while (rays_finished_probe_flag);

      _all_rays_finished = std::accumulate(
          _rays_finished_per_proc.begin(), _rays_finished_per_proc.end(), _rays_finished);

      // Reset
      communication_clicks = 0;
    }

    // Send messages about rays being finished
    for (processor_id_type i = 0; i < num_procs; i++)
    {
      if (i != pid)
      {
        if ((!_rays_finished_requests_sent[i] || _rays_finished_requests[i].test()) &&
            _rays_finished > _rays_finished_requests_temps[i])
        {
          _rays_finished_requests_temps[i] = _rays_finished;
          _comm.send(i,
                     _rays_finished_requests_temps[i],
                     _rays_finished_requests[i],
                     _rays_finished_requests_tag);
          _rays_finished_requests_sent[i] = true;
        }
      }
    }

    if (rays_started_request.test())
      if (_all_rays_started == _all_rays_finished)
        return;

    communication_clicks++;
  }
}

void
RayTracingStudy::bsPropagate()
{
  // std::cout << "Using BS!!!" << std::endl;

  Parallel::Request rays_started_request;
  Parallel::Request rays_finished_request;

  // Get the number of rays that were started in the whole domain
  nonblockingSum(_comm, _rays_started, _all_rays_started, rays_started_request);

  // Keep bouncing the rays around until they've all completed
  while (true)
  {
    bool receiving = false;
    bool sending = false;

    Parallel::Request some_left_request;
    unsigned int some_left = 0;
    unsigned int all_some_left = 1;

    do
    {
      _receive_buffer.receive();
      flushBuffers();

      receiving = _receive_buffer.currentlyReceiving();

      sending = false;
      for (auto & send_buffer : _send_buffers)
        sending = sending || send_buffer.second->currentlySending() ||
                  send_buffer.second->currentlyBuffered();

      if (!receiving && !sending && some_left_request.test() && all_some_left)
      {
        some_left = receiving || sending;

        nonblockingSum(_comm, some_left, all_some_left, some_left_request);
      }
    } while (receiving || sending || !some_left_request.test() || all_some_left);

    _comm.barrier();

    receiveAndTrace();

    if (rays_started_request.test() && rays_finished_request.test())
    {
      if (_all_rays_started == _all_rays_finished)
        return;

      _rays_finished_temp = _rays_finished;
      nonblockingSum(_comm, _rays_finished_temp, _all_rays_finished, rays_finished_request);
    }
  }
}
