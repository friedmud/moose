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

#ifndef ADDNEWSTEPPERACTION_H
#define ADDNEWSTEPPERACTION_H

#include "MooseObjectAction.h"

//Forward Declaration
class AddNewStepperAction;

template<>
InputParameters validParams<AddNewStepperAction>();


class AddNewStepperAction : public MooseObjectAction
{
public:
  AddNewStepperAction(InputParameters params);

  virtual void act() override;
};

#endif // ADDNEWSTEPPERACTION_H
