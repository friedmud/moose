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

#ifndef ADDRayBOUNDARYCONDITIONACTION_H
#define ADDRayBOUNDARYCONDITIONACTION_H

#include "MooseObjectAction.h"

class AddRayBoundaryConditionAction;

template <>
InputParameters validParams<AddRayBoundaryConditionAction>();

class AddRayBoundaryConditionAction : public MooseObjectAction
{
public:
  AddRayBoundaryConditionAction(InputParameters params);

  virtual void act();
};

#endif // ADDRayBOUNDARYCONDITIONACTION_H
