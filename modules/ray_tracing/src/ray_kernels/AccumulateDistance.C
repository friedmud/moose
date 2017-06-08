/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "AccumulateDistance.h"

// Local Includes
#include "RayProblem.h"
#include "RaySystem.h"

// MOOSE Includes
#include "MooseMesh.h"

template <>
InputParameters
validParams<AccumulateDistance>()
{
  InputParameters params = validParams<RayKernel>();

  return params;
}

AccumulateDistance::AccumulateDistance(const InputParameters & params) : RayKernel(params) {}

AccumulateDistance::~AccumulateDistance() {}

void
AccumulateDistance::onSegment(const Elem * /*elem*/,
                              const Point & start,
                              const Point & end,
                              bool ends_in_elem)
{
  _ray_data[0] += (end - start).norm();

  // Let's dump the accumulated distance into a field in the final element
  if (ends_in_elem)
    _group_solution_values[_current_offset] += _ray_data[0];
}

void
AccumulateDistance::setRay(const std::shared_ptr<Ray> & ray)
{
  RayKernel::setRay(ray);

  _ray_data = &_ray->data()[0];
}
