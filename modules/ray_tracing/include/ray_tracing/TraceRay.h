#ifndef TRACERAY_H
#define TRACERAY_H

// MOOSE Includes
#include "MooseTypes.h"

// forward declares
class Ray;
class MooseMesh;
class FEProblem;
class RayKernel;
class RayProblemBase;

namespace TraceRay
{
void traceRay(const std::shared_ptr<Ray> & ray,
              RayProblemBase & ray_problem,
              const MooseMesh & mesh,
              unsigned int halo_size,
              Real ray_max_distance,
              Real ray_length,
              THREAD_ID tid);
}

#endif
