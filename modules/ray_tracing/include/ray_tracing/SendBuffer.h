#ifndef SENDBUFFER_H
#define SENDBUFFER_H

// Local Includes
#include "Ray.h"
#include "RayTracingMethod.h"
#include "RayTracingCommon.h"

// libMesh Includes
#include "libmesh/parallel.h"
#include "libmesh/parallel_object.h"

// System Includes
#include <list>

/**
 * Controls the Ray buffer headed _one_ other processor
 *
 * Will automatically send the buffer when it reaches max size
 */
class SendBuffer : public ParallelObject
{
public:
  /**
   * Constructor
   * @param comm The Communicator to use
   * @param pid The processor ID this buffer will send to
   */
  SendBuffer(const Parallel::Communicator & comm,
             processor_id_type pid,
             unsigned int max_buff_size,
             RayTracingMethod method)
    : ParallelObject(comm), _pid(pid), _max_buff_size(max_buff_size), _method(method)
  {
    _buffer.reserve(_max_buff_size);
  }

  /**
   * Destructor: ensures that all send requests have completed
   */
  ~SendBuffer()
  {
    for (auto & request : _requests)
      request->wait();

    cleanupRequests();
  }

  /**
   * Set the current algorithm
   */
  void setMethod(RayTracingMethod method) { _method = method; }

  /**
   * Whether or not messges are currently being sent
   */
  bool currentlySending() { return _requests.size(); }

  /**
   * Whether or not rays are currently waiting to be sent
   */
  bool currentlyBuffered() { return _buffer.size(); }

  /**
   * Add a Ray to the buffer.  May cause the buffer
   * to be communicated.
   */
  void addRay(std::shared_ptr<Ray> & ray)
  {
    _buffer.push_back(ray);

    // 2000 is a safety net to try to stay below 1MB
    if ((_method == SMART && _buffer.size() == _max_buff_size) || (_buffer.size() > 2000))
      forceSend();
  }

  /**
   * Forces a Send for all currently buffered Rays
   */
  void forceSend()
  {
    if (!_buffer.empty())
    {
      auto req = std::make_shared<Parallel::Request>();

      _requests.push_back(req);

      _communicator.nonblocking_send_packed_range(
          _pid, &_buffer, _buffer.begin(), _buffer.end(), *req, _ray_buffer_tag);

      _buffer.clear();
      _buffer.reserve(_max_buff_size);
    }

    cleanupRequests();
  }

  /**
   * Wait for all requests to finish
   */
  void waitAll()
  {
    for (auto & request : _requests)
      request->wait();
  }

  /**
   * Clear all existing data
   */
  void clear()
  {
    /// The buffer
    _buffer.clear();

    /// List of Requests
    _requests.clear();
  }

  /**
   * Checks to see if any Requests can be finished
   */
  void cleanupRequests()
  {
    _requests.remove_if([](std::shared_ptr<Parallel::Request> & req) {
      if (req->test())
      {
        req->wait(); // MUST call wait() to do post_wait_work
        return true;
      }
      else
        return false;
    });
  }

protected:
  /// Current size of the send buffer (double on purpose!)
  double _current_size;

  /// The processor ID this SendBuffer will send to
  processor_id_type _pid;

  /// Maximum size of the buffer (in Rays)
  unsigned int _max_buff_size;

  /// The buffer
  std::vector<std::shared_ptr<Ray>> _buffer;

  /// List of Requests
  std::list<std::shared_ptr<Parallel::Request>> _requests;

  /// MessageTag for sending rays
  Parallel::MessageTag _ray_buffer_tag = Parallel::MessageTag(RAY_BUFFER);

  /// The current algorithm
  RayTracingMethod _method;
};

#endif
