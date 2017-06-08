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

#ifndef ADDRayMATERIALACTION_H
#define ADDRayMATERIALACTION_H

#include "MooseObjectAction.h"

class AddRayMaterialAction;

template <>
InputParameters validParams<AddRayMaterialAction>();

class AddRayMaterialAction : public MooseObjectAction
{
public:
  AddRayMaterialAction(InputParameters params);

  virtual void act();
};

#endif // ADDRayMATERIALACTION_H
