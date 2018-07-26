#ifndef RAYTRACINGMETHOD_H
#define RAYTRACINGMETHOD_H

#include "LIFOBuffer.h"
#include "CircularBuffer.h"

class Ray;

enum RayTracingMethod
{
  SMART,
  HARM,
  BS
};

/**
 * Type of buffer being used for the work buffer
 */
typedef MooseUtils::CircularBuffer<std::shared_ptr<Ray>> WorkBufferType;
// typedef MooseUtils::LIFOBuffer<std::shared_ptr<Ray>> WorkBufferType;

#endif
