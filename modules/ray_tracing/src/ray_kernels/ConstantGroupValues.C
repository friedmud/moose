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

#include "ConstantGroupValues.h"

// Local Includes
#include "RayProblem.h"
#include "RaySystem.h"

// MOOSE Includes
#include "MooseMesh.h"

template <>
InputParameters
validParams<ConstantGroupValues>()
{
  InputParameters params = validParams<RayKernel>();

  params.addRequiredParam<std::vector<Real>>("group_values", "The value for each group");

  return params;
}

ConstantGroupValues::ConstantGroupValues(const InputParameters & params)
  : RayKernel(params),
    _num_groups(_ray_problem.numGroups()),
    _input_group_values(getParam<std::vector<Real>>("group_values"))
{
  if (_num_groups != _input_group_values.size())
    mooseError("group_values size does not match num_groups! In ", name());
}

ConstantGroupValues::~ConstantGroupValues() {}

void
ConstantGroupValues::onSegment(const Elem * /*elem*/,
                               const Point & /* start */,
                               const Point & /* end */,
                               bool /* ends_in_elem */)
{
  for (unsigned int g = 0; g < _num_groups; g++)
    _group_solution_values[_current_offset + g] = _input_group_values[g];
}

void
ConstantGroupValues::setRay(const std::shared_ptr<Ray> & ray)
{
  RayKernel::setRay(ray);
}
