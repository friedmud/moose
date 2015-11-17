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

#include "ExecuteInitialExactGeometry.h"
#include "FEProblem.h"

template<>
InputParameters validParams<ExecuteInitialExactGeometry>()
{
  InputParameters params = validParams<Action>();
  return params;
}


ExecuteInitialExactGeometry::ExecuteInitialExactGeometry(InputParameters params) :
    Action(params)
{
}

void
ExecuteInitialExactGeometry::act()
{
  _problem->executeInitialExactGeometry();
}
