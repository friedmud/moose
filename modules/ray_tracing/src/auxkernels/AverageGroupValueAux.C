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

#include "AverageGroupValueAux.h"

template <>
InputParameters
validParams<AverageGroupValueAux>()
{
  InputParameters params = validParams<GroupValueAuxBase>();
  return params;
}

AverageGroupValueAux::AverageGroupValueAux(const InputParameters & parameters)
  : GroupValueAuxBase(parameters)
{
}

Real
AverageGroupValueAux::computeValue()
{
  Real total = 0;

  for (unsigned int g = 0; g < _num_groups; g++)
    total += _group_values[_current_offset + g];

  return total / static_cast<Real>(_num_groups);
}
