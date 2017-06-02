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

#include "RayBC.h"

// Local Includes
#include "RunStudy.h"
#include "RAYProblem.h"
#include "RAYSystem.h"

// MOOSE Includes
#include "MooseMesh.h"

template <>
InputParameters
validParams<RayBC>()
{
  InputParameters params = validParams<MooseObject>();

  params += validParams<SetupInterface>();
  params += validParams<BoundaryRestrictableRequired>();

  params.registerBase("RayBC");

  return params;
}

RayBC::RayBC(const InputParameters & params)
  : MooseObject(params),
    SetupInterface(this),
    Restartable(params, "RayBCs"),
    BoundaryRestrictableRequired(params, false), // false for sidesets
    Coupleable(this, false),
    GeometricSearchInterface(this),
    _ray_problem(dynamic_cast<RAYProblem &>(*params.get<FEProblem *>("_fe_problem"))),
    _ray_sys(_ray_problem.raySystem()),
    _tid(params.get<THREAD_ID>("_tid"))
{
}

RayBC::~RayBC() {}
