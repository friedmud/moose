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

#include "RayMaterial.h"

// Local Includes
#include "RayProblem.h"

// MOOSE Includes
#include "MooseMesh.h"

template <>
InputParameters
validParams<RayMaterial>()
{
  InputParameters params = validParams<MooseObject>();

  params += validParams<SetupInterface>();
  params += validParams<BlockRestrictable>();

  params.addParam<bool>("fissionable", false, "Whether or not this material produces fissions");

  params.registerBase("RayMaterial");

  return params;
}

RayMaterial::RayMaterial(const InputParameters & params)
  : MooseObject(params),
    SetupInterface(this),
    Restartable(params, "RayMaterials"),
    BlockRestrictable(params),
    Coupleable(this, false),
    UserObjectInterface(this),
    TransientInterface(this),
    _ray_problem(dynamic_cast<RayProblemBase &>(*params.get<FEProblem *>("_fe_problem"))),
    _num_groups(_ray_problem.numGroups()),
    _tid(params.get<THREAD_ID>("_tid"))
{
}
