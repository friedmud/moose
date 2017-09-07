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

#include "AddRayBCAction.h"

// Local Includes
#include "RayProblem.h"

template <>
InputParameters
validParams<AddRayBoundaryConditionAction>()
{
  return validParams<MooseObjectAction>();
}

AddRayBoundaryConditionAction::AddRayBoundaryConditionAction(InputParameters params)
  : MooseObjectAction(params)
{
}

void
AddRayBoundaryConditionAction::act()
{
  dynamic_cast<RayProblemBase *>(_problem.get())
      ->raySystem()
      .addRayBC(_type, _name, _moose_object_pars);
}
