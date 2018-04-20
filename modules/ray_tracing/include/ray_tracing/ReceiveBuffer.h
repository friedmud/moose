#ifndef RECEIVEBUFFER_H
#define RECEIVEBUFFER_H

// Local Includes
#include "Ray.h"
#include "RayTracingMethod.h"

// Moose Includes
#include "MooseError.h"
#include "CircularBuffer.h"

// libMesh Includes
#include "libmesh/parallel.h"
#include "libmesh/parallel_object.h"
#include "libmesh/mesh.h"

// System Includes
#include <list>
#include <queue>

/**
 * Keeps track of many non-blocking receives
 */
class ReceiveBuffer : public ParallelObject
{
public:
  /**
   * Constructor
   * @param comm The Communicator to use
   */
  ReceiveBuffer(const Parallel::Communicator & comm,
                MeshBase & mesh,
                unsigned int clicks_per_receive,
                bool blocking = false)
    : ParallelObject(comm),
      _mesh(mesh),
      _clicks_per_receive(clicks_per_receive),
      _blocking(blocking),
      _method(SMART)
  {
  }

  /**
   * Destructor: ensures that all send requests have completed
   */
  ~ReceiveBuffer()
  {
    if (!_requests.empty())
      mooseError("Some requests not serviced!");
  }

  /**
   * Set the current algorithm
   */
  void setMethod(RayTracingMethod method) { _method = method; }

  /**
   * Whether or not there are messages that are being currently received
   */
  bool currentlyReceiving() { return _requests.size(); }

  /**
   * The number of rays received since the last reset
   */
  unsigned long int raysReceived() { return _rays_received; }

  /**
   * The number of buffers received since the last reset
   */
  unsigned long int buffersReceived() { return _buffers_received; }

  /**
   * Start receives for all currently available messages
   *
   * Adds the to the working buffer
   */
  void receive(MooseUtils::CircularBuffer<std::shared_ptr<Ray>> & working_buffer)
  {
    bool flag = false;
    Parallel::Status stat;

    static unsigned int current_clicks = 0;

    if (current_clicks % _clicks_per_receive == 0)
    {
      current_clicks = 0;

      // Receive and process a bunch of Rays
      do
      {
        stat = _communicator.template packed_range_probe<std::shared_ptr<Ray>>(
            Parallel::any_source, _ray_buffer_tag, flag);

        if (flag)
        {
          auto req = std::make_shared<Parallel::Request>();
          auto rays = std::make_shared<std::vector<std::shared_ptr<Ray>>>();

          if (_method == HARM || _method == BS)
            blocking_receive_packed_range(_communicator,
                                          stat.source(),
                                          &_mesh,
                                          std::back_inserter(*rays),
                                          (std::shared_ptr<Ray> *)(libmesh_nullptr),
                                          *req,
                                          stat,
                                          _ray_buffer_tag);
          else
            _communicator.nonblocking_receive_packed_range(
                stat.source(),
                &_mesh,
                std::back_inserter(*rays),
                (std::shared_ptr<Ray> *)(libmesh_nullptr),
                *req,
                stat,
                _ray_buffer_tag);

          _requests.emplace_back(std::make_pair(req, rays));
        }
      } while (flag);
    }

    current_clicks++;

    cleanupRequests(working_buffer);
  }

  /**
   * Wait for all requests to finish
   */
  void waitAll(MooseUtils::CircularBuffer<std::shared_ptr<Ray>> & working_buffer)
  {
    for (auto & request_pair : _requests)
    {
      request_pair.first->wait();
      working_buffer.append(request_pair.second->begin(), request_pair.second->end());
    }
  }

  /**
   * Clear all existing data
   */
  void clear()
  {
    // List of Requests and buffers for each request
    _requests.clear();

    _receive_loop_time = std::chrono::steady_clock::duration::zero();
    _cleanup_requests_time = std::chrono::steady_clock::duration::zero();

    _rays_received = 0;
    _buffers_received = 0;
  }

  /**
   * Number of iterations before looking for new rays
   */
  unsigned int receiveClicks() { return _clicks_per_receive; }

  /**
   * Set Number of iterations before looking for new rays
   */
  void setReceiveClicks(unsigned int clicks_per_receive)
  {
    _clicks_per_receive = clicks_per_receive;
  }

  /**
   * Amount of time (in seconds) spent in the loop checking for messages and creating Receives
   */
  Real receiveLoopTime() { return std::chrono::duration<Real>(_receive_loop_time).count(); }

  /**
   * Amount of time (in seconds) spent finishing receives and collecting rays
   */
  Real cleanupRequestsTime() { return std::chrono::duration<Real>(_cleanup_requests_time).count(); }

  /**
   * Checks to see if any Requests can be finished
   */
  void cleanupRequests(MooseUtils::CircularBuffer<std::shared_ptr<Ray>> & working_buffer)
  {
    //    auto cleanup_requests_start = std::chrono::steady_clock::now();

    _requests.remove_if([&](std::pair<std::shared_ptr<Parallel::Request>,
                                      std::shared_ptr<std::vector<std::shared_ptr<Ray>>>> &
                                request_pair) {
      auto req = request_pair.first;
      auto rays = request_pair.second;

      if (req->test()) // See if the receive has completed
      {
        req->wait(); // MUST call wait() to do post_wait_work which actually fills the ray buffer

        _buffers_received++;
        _rays_received += rays->size();

        working_buffer.append(rays->begin(), rays->end());

        return true;
      }
      else
        return false;
    });

    //    _cleanup_requests_time += std::chrono::steady_clock::now() - cleanup_requests_start;
  }

protected:
  /// The Mesh we're working on
  MeshBase & _mesh;

  /// Number of iterations to wait before looking for more rays
  unsigned int _clicks_per_receive;

  /// List of Requests and buffers for each request
  std::list<std::pair<std::shared_ptr<Parallel::Request>,
                      std::shared_ptr<std::vector<std::shared_ptr<Ray>>>>>
      _requests;

  /// MessageTag for sending rays
  Parallel::MessageTag _ray_buffer_tag = Parallel::MessageTag(RAY_BUFFER);

  /// Wether or not to use blocking receives
  bool _blocking;

  /// Receive loop time
  std::chrono::steady_clock::duration _receive_loop_time;

  /// Time cleaning up requests
  std::chrono::steady_clock::duration _cleanup_requests_time;

  /// The current algorithm
  RayTracingMethod _method;

  /// Total rays received
  unsigned long int _rays_received;

  /// Total ray buffers received
  unsigned long int _buffers_received;

private:
  template <typename Context, typename OutputIter, typename T>
  inline void blocking_receive_packed_range(const Parallel::Communicator & comm,
                                            const unsigned int src_processor_id,
                                            Context * context,
                                            OutputIter out,
                                            const T * /* output_type */,
                                            Parallel::Request & req,
                                            Parallel::Status & stat,
                                            const Parallel::MessageTag & tag) const
  {
    libmesh_experimental();

    typedef typename Parallel::Packing<T>::buffer_type buffer_t;

    // Receive serialized variable size objects as a sequence of
    // buffer_t.
    // Allocate a buffer on the heap so we don't have to free it until
    // after the Request::wait()
    std::vector<buffer_t> * buffer = new std::vector<buffer_t>(stat.size());
    comm.receive(src_processor_id, *buffer, tag);

    // Make the Request::wait() handle unpacking the buffer
    req.add_post_wait_work(
        new libMesh::Parallel::PostWaitUnpackBuffer<std::vector<buffer_t>, Context, OutputIter, T>(
            *buffer, context, out));

    // Make the Request::wait() then handle deleting the buffer
    req.add_post_wait_work(
        new libMesh::Parallel::PostWaitDeleteBuffer<std::vector<buffer_t>>(buffer));
  }
};

#endif
