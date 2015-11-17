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

#ifndef ADDEXACTGEOMETRYACTION_H
#define ADDEXACTGEOMETRYACTION_H

#include "MooseObjectAction.h"

class AddExactGeometryAction;

template<>
InputParameters validParams<AddExactGeometryAction>();


class AddExactGeometryAction : public MooseObjectAction
{
public:
  AddExactGeometryAction(InputParameters params);

  virtual void act();
};

#endif // ADDEXACTGEOMETRYACTION_H
