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

#ifndef EXECUTEINITIALEXACTGEOMETRY_H
#define EXECUTEINITIALEXACTGEOMETRY_H

#include "Action.h"

class ExecuteInitialExactGeometry;

template<>
InputParameters validParams<ExecuteInitialExactGeometry>();


class ExecuteInitialExactGeometry : public Action
{
public:
  ExecuteInitialExactGeometry(InputParameters params);

  virtual void act();
};

#endif // EXECUTEINITIALEXACTGEOMETRY_H
