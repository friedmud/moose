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

#include "GroupValueAuxBase.h"

template <>
InputParameters
validParams<GroupValueAuxBase>()
{
  InputParameters params = validParams<AuxKernel>();
  return params;
}

GroupValueAuxBase::GroupValueAuxBase(const InputParameters & parameters)
  : AuxKernel(parameters),
    _ray_problem(dynamic_cast<RayProblemBase &>(*parameters.get<FEProblem *>("_fe_problem"))),
    _ray_sys(_ray_problem.raySystem()),
    _num_groups(_ray_problem.numGroups()),
    _current_offset(_ray_sys.currentOffset(_tid)),
    _group_values(_ray_sys.currentGroupSolutionValues())
{
}
