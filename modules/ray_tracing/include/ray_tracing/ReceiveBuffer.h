#ifndef RECEIVEBUFFER_H
#define RECEIVEBUFFER_H

// Local Includes
#include "Ray.h"
#include "RayTracingMethod.h"

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
  ~ReceiveBuffer() { waitAll(); }

  /**
   * Set the current algorithm
   */
  void setMethod(RayTracingMethod method) { _method = method; }

  /**
   * Whether or not there are messages that are being currently received
   */
  bool currentlyReceiving() { return _requests.size(); }

  /**
   * Start receives for all currently available messages
   */
  void receive()
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

          _requests.push_back(std::make_pair(req, rays));
        }
      } while (flag);
    }

    current_clicks++;

    cleanupRequests();
  }

  /**
   * Grab a vector of Rays to work on
   *
   * @return A vector of Rays or NULL if there are none
   */
  std::shared_ptr<std::vector<std::shared_ptr<Ray>>> getRays()
  {
    if (_available_rays.empty())
      return NULL;

    unsigned int total_size = 0;

    // Get the total count of Rays
    for (auto & ray_vec : _available_rays)
      total_size += ray_vec->size();

    auto rays = std::make_shared<std::vector<std::shared_ptr<Ray>>>();

    rays->reserve(total_size);

    for (auto & ray_vec : _available_rays)
      rays->insert(rays->end(), ray_vec->begin(), ray_vec->end());

    _available_rays.clear();

    return rays;
  }

  /**
   * Wait for all requests to finish
   */
  void waitAll()
  {
    for (auto & request_pair : _requests)
    {
      request_pair.first->wait();
      _available_rays.push_back(request_pair.second);
    }
  }

  /**
   * Clear all existing data
   */
  void clear()
  {
    // List of Requests and buffers for each request
    _requests.clear();

    // Available buffers of Rays that have been received
    // Use the "swap trick"
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<Ray>>>> empty;
    std::swap(_available_rays, empty);

    _receive_loop_time = std::chrono::steady_clock::duration::zero();
    _cleanup_requests_time = std::chrono::steady_clock::duration::zero();
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

protected:
  /**
   * Checks to see if any Requests can be finished
   */
  void cleanupRequests()
  {
    //    auto cleanup_requests_start = std::chrono::steady_clock::now();

    _requests.remove_if([&](
        std::pair<std::shared_ptr<Parallel::Request>,
                  std::shared_ptr<std::vector<std::shared_ptr<Ray>>>> & request_pair) {
      auto req = request_pair.first;
      auto rays = request_pair.second;

      if (req->test()) // See if the receive has completed
      {
        req->wait(); // MUST call wait() to do post_wait_work which actually fills the ray buffer

        _available_rays.push_back(rays);

        return true;
      }
      else
        return false;
    });

    //    _cleanup_requests_time += std::chrono::steady_clock::now() - cleanup_requests_start;
  }

  /// The Mesh we're working on
  MeshBase & _mesh;

  /// Number of iterations to wait before looking for more rays
  unsigned int _clicks_per_receive;

  /// List of Requests and buffers for each request
  std::list<std::pair<std::shared_ptr<Parallel::Request>,
                      std::shared_ptr<std::vector<std::shared_ptr<Ray>>>>>
      _requests;

  /// Available buffers of Rays that have been received
  std::vector<std::shared_ptr<std::vector<std::shared_ptr<Ray>>>> _available_rays;

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
