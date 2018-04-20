#ifndef CIRCULARBUFFER_H
#define CIRCULARBUFFER_H

#include <vector>

namespace MooseUtils
{

/**
 * An optimized circular buffer.
 *
 * This also always ensures that begin() < end().
 * That means that if the end of the buffer capacity is reached an O(N) operation
 * will be used to move data from the end of the capacity to the beginning
 *
 * This is done so that begin()/end() iterators can be used in OpenMP loops
 * Also means that operator[] works sequentially between 0 and size()
 *
 * It will also automatically grow larger if capacity is reached
 */
template <typename T>
class CircularBuffer
{
public:
  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

  /**
   * Create an empty circular buffer
   */
  CircularBuffer();

  /**
   * Create a buffer with a specific capacity
   */
  CircularBuffer(std::size_t capacity);

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
   * Remove the elements in the range [begin(); last)
   *
   * Note that erased items are not guaranteed to be freed immediately
   */
  typename CircularBuffer<T>::iterator erase(iterator last);

  /**
   * Remove the first num elements
   *
   * Note that erased items are not guaranteed to be freed immediately
   */
  typename CircularBuffer<T>::iterator erase(std::size_t num);

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
   * Access the raw underlying array
   *
   * Caution!  Use with care!
   * The values within the underlying array
   * can move and be removed
   */
  std::vector<T> & data() { return _data; }

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
CircularBuffer<T>::CircularBuffer() : _begin_pos(0), _end_pos(0)
{
}

template <typename T>
CircularBuffer<T>::CircularBuffer(std::size_t capacity)
  : _data(capacity), _begin_pos(0), _end_pos(0)
{
}

template <typename T>
void
CircularBuffer<T>::setCapacity(std::size_t capacity)
{
  _data.resize(capacity);
}

template <typename T>
std::size_t
CircularBuffer<T>::capacity()
{
  return _data.size();
}

template <typename T>
void
CircularBuffer<T>::setSize(std::size_t size)
{
  _end_pos = newEnd(_begin_pos + size);
}

template <typename T>
std::size_t
CircularBuffer<T>::size()
{
  return _end_pos - _begin_pos;
}

template <typename T>
void
CircularBuffer<T>::push_back(const T & value)
{
  _end_pos = newEnd(_end_pos + 1);

  _data[_end_pos - 1] = value;
}

template <typename T>
void
CircularBuffer<T>::append(const_iterator in_begin, const_iterator in_end)
{
  auto additional_size = std::distance(in_begin, in_end);

  auto new_end = newEnd(_end_pos + additional_size);

  std::copy(in_begin, in_end, end());

  _end_pos = new_end;
}

template <typename T>
void
CircularBuffer<T>::clear()
{
  _begin_pos = 0;
  _end_pos = 0;
}

template <typename T>
typename CircularBuffer<T>::iterator
CircularBuffer<T>::erase(iterator last)
{
  mooseAssert(last <= end(), "Cannot erase past the last entry!");

  _begin_pos += std::distance(begin(), last);

  // If there's nothing in the buffer - let's reset the positions
  if (_begin_pos == _end_pos)
    clear();

  return begin();
}

template <typename T>
typename CircularBuffer<T>::iterator
CircularBuffer<T>::erase(std::size_t num)
{
  _begin_pos += num;
  mooseAssert(_begin_pos <= _end_pos, "Cannot erase past the last entry!");

  // If there's nothing in the buffer - let's reset the positions
  if (_begin_pos == _end_pos)
    clear();

  return begin();
}

template <typename T>
typename CircularBuffer<T>::iterator
CircularBuffer<T>::begin()
{
  return _data.begin() + _begin_pos;
}

template <typename T>
typename CircularBuffer<T>::const_iterator
CircularBuffer<T>::begin() const
{
  return _data.begin() + _begin_pos;
}

template <typename T>
typename CircularBuffer<T>::iterator
CircularBuffer<T>::end()
{
  return _data.begin() + _end_pos;
}

template <typename T>
typename CircularBuffer<T>::const_iterator
CircularBuffer<T>::end() const
{
  return _data.begin() + _end_pos;
}

template <typename T>
T & CircularBuffer<T>::operator[](std::size_t index)
{
  mooseAssert(_begin_pos + index < _end_pos, "Attempt to access off end of CircularBuffer!");
  return _data[_begin_pos + index];
}

template <typename T>
const T & CircularBuffer<T>::operator[](std::size_t index) const
{
  mooseAssert(_begin_pos + index < _end_pos, "Attempt to access off end of CircularBuffer!");
  return _data[_begin_pos + index];
}

template <typename T>
std::size_t
CircularBuffer<T>::newEnd(std::size_t new_end)
{
  if (new_end > _data.size())
  {
    auto new_size = new_end - _begin_pos;

    // See if we need to grow our capacity
    if (_begin_pos == 0) // If we're already using the beginning - just resize
      _data.resize(new_size * 2);
    else
    {
      // Try to move everything to the beginning of the array
      std::copy(begin(), end(), _data.begin());

      _begin_pos = 0;
      new_end = new_size;

      // If there still isn't room... add space
      if (new_end > _data.size())
        _data.resize(new_size * 2);
    }
  }

  return new_end;
}
}

#endif
