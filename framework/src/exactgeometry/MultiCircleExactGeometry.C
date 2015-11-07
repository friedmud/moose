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

#include "MultiCircleExactGeometry.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "Assembly.h"
#include "MooseVariable.h"

template<>
InputParameters validParams<MultiCircleExactGeometry>()
{
  InputParameters params = validParams<ExactGeometry>();
  return params;
}

MultiCircleExactGeometry::MultiCircleExactGeometry(const InputParameters & parameters) :
    ExactGeometry(parameters)
{
}

void
MultiCircleExactGeometry::moveNode(boundary_id_type boundary_id, Node & node)
{
}
