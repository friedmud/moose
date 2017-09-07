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

#include "AddRayMaterialAction.h"

// Local Includes
#include "RayProblem.h"

template <>
InputParameters
validParams<AddRayMaterialAction>()
{
  return validParams<MooseObjectAction>();
}

AddRayMaterialAction::AddRayMaterialAction(InputParameters params) : MooseObjectAction(params) {}

void
AddRayMaterialAction::act()
{
  dynamic_cast<RayProblemBase *>(_problem.get())
      ->raySystem()
      .addRayMaterial(_type, _name, _moose_object_pars);
}
