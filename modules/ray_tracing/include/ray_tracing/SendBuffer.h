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
#include <algorithm>

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
             unsigned int min_buff_size,
             Real buffer_growth_multiplier,
             Real buffer_shrink_multiplier,
             RayTracingMethod method)
    : ParallelObject(comm),
      _pid(pid),
      _max_buff_size(max_buff_size),
      _min_buff_size(min_buff_size),
      _buffer_growth_multiplier(buffer_growth_multiplier),
      _buffer_shrink_multiplier(buffer_shrink_multiplier),
      _current_buff_size(_min_buff_size),
      _current_buff_size_real(_min_buff_size),
      _method(method),
      _rays_sent(0),
      _buffers_sent(0)
  {
    std::cout << "Min buff size: " << _min_buff_size << std::endl;

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
   * Get the number of rays sent from this buffer
   */
  unsigned long int raysSent() { return _rays_sent; }

  /**
   * Get the number of buffers sent from this buffer
   */
  unsigned long int buffersSent() { return _buffers_sent; }

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
    if ((_method == SMART &&
         (_buffer.size() >= _current_buff_size || _buffer.size() == _max_buff_size)) ||
        (_buffer.size() > 200000))
    {
      //      std::cout << "Sending " << _buffer.size() << " rays to " << _pid << "\n";

      _current_buff_size_real = std::min(_buffer_growth_multiplier * _current_buff_size_real,
                                         static_cast<Real>(_max_buff_size));

      if (_current_buff_size != static_cast<unsigned int>(_current_buff_size_real))
      {
        _current_buff_size = static_cast<unsigned int>(_current_buff_size_real);
        std::cout << "Increasing buffer size to: " << _current_buff_size << "\n";
      }

      forceSend(false);
    }
  }

  /**
   * Forces a Send for all currently buffered Rays
   */
  void forceSend(bool shrink_current_buff_size = true)
  {
    if (!_buffer.empty())
    {
      auto req = std::make_shared<Parallel::Request>();

      _requests.push_back(req);

      _rays_sent += _buffer.size();
      _buffers_sent++;

      _communicator.nonblocking_send_packed_range(
          _pid, &_buffer, _buffer.begin(), _buffer.end(), *req, _ray_buffer_tag);

      _buffer.clear();
      _buffer.reserve(_max_buff_size);

      if (_method == SMART && shrink_current_buff_size)
      {
        _current_buff_size_real = std::max(static_cast<Real>(_min_buff_size),
                                           _current_buff_size_real * _buffer_shrink_multiplier);

        if (_current_buff_size != static_cast<unsigned int>(_current_buff_size_real))
        {
          _current_buff_size = static_cast<unsigned int>(_current_buff_size_real);
          std::cout << "Cutting buffer to: " << _current_buff_size << "\n";
        }
      }
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

    // Reset Counters
    _rays_sent = 0;
    _buffers_sent = 0;
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

  /// Minimum size of the buffer (in Rays)
  unsigned int _min_buff_size;

  /// Multiplier for the buffer size for growing the buffer
  Real _buffer_growth_multiplier;

  /// Multiplier for the buffer size for shrinking the buffer
  Real _buffer_shrink_multiplier;

  /// Current size of the buffer (in Rays)
  unsigned int _current_buff_size;

  /// Running buffer size
  Real _current_buff_size_real;

  /// The buffer
  std::vector<std::shared_ptr<Ray>> _buffer;

  /// List of Requests
  std::list<std::shared_ptr<Parallel::Request>> _requests;

  /// MessageTag for sending rays
  Parallel::MessageTag _ray_buffer_tag = Parallel::MessageTag(RAY_BUFFER);

  /// The current algorithm
  RayTracingMethod _method;

  /// Counters
  unsigned long int _rays_sent;
  unsigned long int _buffers_sent;
};

#endif
