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

#ifndef ADDRayKERNELACTION_H
#define ADDRayKERNELACTION_H

#include "MooseObjectAction.h"

class AddRayKernelAction;

template <>
InputParameters validParams<AddRayKernelAction>();

class AddRayKernelAction : public MooseObjectAction
{
public:
  AddRayKernelAction(InputParameters params);

  virtual void act();
};

#endif // ADDRayKERNELACTION_H
