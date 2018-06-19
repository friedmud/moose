#ifndef LIFOBUFFER_H
#define LIFOBUFFER_H

#include "MooseError.h"

#include <vector>

namespace MooseUtils
{

/**
 * An optimized LIFO (Last In First Out) buffer.
 *
 * Begin()/end() iterators can be used in OpenMP loops
 * Also operator[] works sequentially between 0 and size()
 *
 * It will also automatically grow larger if capacity is reached
 */
template <typename T>
class LIFOBuffer
{
public:
  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

  /**
   * Create an empty circular buffer
   */
  LIFOBuffer();

  /**
   * Create a buffer with a specific capacity
   */
  LIFOBuffer(std::size_t capacity);

  /**
   * Resize the capacity
   */
  void setCapacity(std::size_t capacity);

  /**
   * Get the capacity
   */
  std::size_t capacity();

  /**
   * Set the size
   */
  void setSize(std::size_t size);

  /**
   * Get the size
   */
  std::size_t size();

  /**
   * Empty?
   */
  bool empty() { return !size(); }

  /**
   * Add a new entry on the end
   */
  void push_back(const T & value);

  /**
   * Add a new entries to the end
   */
  void append(const std::vector<T> & vals) { append(vals.begin(), vals.end()); }

  /**
   * Add a new entries to the end
   *
   * Everything in [in_begin, in_end) is appended
   */
  void append(const_iterator in_begin, const_iterator in_end);

  /**
   * Remove all entries (does not change the capacity)
   * Note: this does NOT at all free any entries
   */
  void clear();

  /**
   * Remove the elements in the range [begin; end())
   *
   * Note that erased items are not guaranteed to be freed immediately
   */
  typename LIFOBuffer<T>::iterator erase(iterator begin);

  /**
   * Remove the last num elements
   *
   * Note that erased items are not guaranteed to be freed immediately
   */
  typename LIFOBuffer<T>::iterator erase(std::size_t num);

  /**
   * Iterator for the first entry in the buffer
   */
  iterator begin();
  const_iterator begin() const;

  /**
   * Iterator for the last entry
   */
  iterator end();
  const_iterator end() const;

  /**
   * Access an entry at index
   */
  T & operator[](std::size_t index);
  const T & operator[](std::size_t index) const;

  /**
   * Use in_data as our data vector
   */
  void swap(std::vector<T> & in_data)
  {
    std::swap(in_data, _data);
    _begin_pos = 0;
    _end_pos = _data.size();
  }

  /**
   * Access the raw underlying array
   *
   * Caution!  Use with care!
   * The values within the underlying array
   * can move and be removed
   */
  std::vector<T> & data() { return _data; }

  /**
   * Expert interface: the current beginning position
   */
  size_t beginPos() { return _begin_pos; }

  /**
   * Expert interface: the current ending position
   */
  size_t endPos() { return _end_pos; }

protected:
  /**
   * Deal with new end
   *
   * @param new_end the proposed new_end position
   * @return The actual position of the new ending
   */
  std::size_t newEnd(std::size_t new_end);

  std::vector<T> _data;

  /// Save the beginning and ending as positions internally
  /// If these were saved as iterators then they would be invalidated
  /// when the storage was resized
  size_t _begin_pos;
  size_t _end_pos;
};

template <typename T>
LIFOBuffer<T>::LIFOBuffer() : _begin_pos(0), _end_pos(0)
{
}

template <typename T>
LIFOBuffer<T>::LIFOBuffer(std::size_t capacity) : _data(capacity), _begin_pos(0), _end_pos(0)
{
}

template <typename T>
void
LIFOBuffer<T>::setCapacity(std::size_t capacity)
{
  _data.resize(capacity);
}

template <typename T>
std::size_t
LIFOBuffer<T>::capacity()
{
  return _data.size();
}

template <typename T>
void
LIFOBuffer<T>::setSize(std::size_t size)
{
  _end_pos = newEnd(_begin_pos + size);
}

template <typename T>
std::size_t
LIFOBuffer<T>::size()
{
  return _end_pos - _begin_pos;
}

template <typename T>
void
LIFOBuffer<T>::push_back(const T & value)
{
  _end_pos = newEnd(_end_pos + 1);

  _data[_end_pos - 1] = value;
}

template <typename T>
void
LIFOBuffer<T>::append(const_iterator in_begin, const_iterator in_end)
{
  auto additional_size = std::distance(in_begin, in_end);

  _end_pos = newEnd(_end_pos + additional_size);

  std::copy(in_begin, in_end, end() - additional_size);
}

template <typename T>
void
LIFOBuffer<T>::clear()
{
  _begin_pos = 0;
  _end_pos = 0;
}

template <typename T>
typename LIFOBuffer<T>::iterator
LIFOBuffer<T>::erase(iterator begin)
{
  mooseAssert(begin <= end(), "Cannot erase past the last entry!");

  _end_pos -= std::distance(begin, end());

  // If there's nothing in the buffer - let's reset the positions
  if (_begin_pos == _end_pos)
    clear();

  return end();
}

template <typename T>
typename LIFOBuffer<T>::iterator
LIFOBuffer<T>::erase(std::size_t num)
{
  mooseAssert(num <= _end_pos, "Cannot erase past the last entry!");

  _end_pos -= num;

  // If there's nothing in the buffer - let's reset the positions
  if (_begin_pos == _end_pos)
    clear();

  return end();
}

template <typename T>
typename LIFOBuffer<T>::iterator
LIFOBuffer<T>::begin()
{
  return _data.begin() + _begin_pos;
}

template <typename T>
typename LIFOBuffer<T>::const_iterator
LIFOBuffer<T>::begin() const
{
  return _data.begin() + _begin_pos;
}

template <typename T>
typename LIFOBuffer<T>::iterator
LIFOBuffer<T>::end()
{
  return _data.begin() + _end_pos;
}

template <typename T>
typename LIFOBuffer<T>::const_iterator
LIFOBuffer<T>::end() const
{
  return _data.begin() + _end_pos;
}

template <typename T>
T & LIFOBuffer<T>::operator[](std::size_t index)
{
  mooseAssert(_begin_pos + index < _end_pos, "Attempt to access off end of LIFOBuffer!");
  return _data[_begin_pos + index];
}

template <typename T>
const T & LIFOBuffer<T>::operator[](std::size_t index) const
{
  mooseAssert(_begin_pos + index < _end_pos, "Attempt to access off end of LIFOBuffer!");
  return _data[_begin_pos + index];
}

template <typename T>
std::size_t
LIFOBuffer<T>::newEnd(std::size_t new_end)
{
  if (new_end > _data.size())
  {
    auto new_size = new_end - _begin_pos;

    _data.resize(std::max(new_size * 2, _data.size() * 2));
  }

  return new_end;
}
}

#endif
