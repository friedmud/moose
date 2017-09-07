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

#include "RayKernel.h"

// Local Includes
#include "RayProblem.h"
#include "RaySystem.h"

// MOOSE Includes
#include "MooseMesh.h"

template <>
InputParameters
validParams<RayKernel>()
{
  InputParameters params = validParams<MooseObject>();

  params += validParams<SetupInterface>();
  params += validParams<BlockRestrictable>();

  params.registerBase("RayKernel");

  return params;
}

RayKernel::RayKernel(const InputParameters & params)
  : MooseObject(params),
    SetupInterface(this),
    Restartable(params, "RayKernels"),
    BlockRestrictable(params),
    Coupleable(this, false),
    _ray_problem(dynamic_cast<RayProblemBase &>(*params.get<FEProblem *>("_fe_problem"))),
    _ray_sys(_ray_problem.raySystem()),
    _tid(params.get<THREAD_ID>("_tid")),
    _num_groups(_ray_problem.numGroups()),
    _current_offset(_ray_sys.currentOffset(_tid)),
    _group_solution_values(_ray_sys.groupSolutionValues(_tid))
{
}

RayKernel::~RayKernel() {}
