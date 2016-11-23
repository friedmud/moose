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

#ifndef CONSTANTNEWSTEPPER_H
#define CONSTANTNEWSTEPPER_H

#include "NewStepper.h"

class ConstantNewStepper;

template<>
InputParameters validParams<ConstantNewStepper>();

class ConstantNewStepper : public NewStepper
{
public:
  ConstantNewStepper(const InputParameters & parameters);
  virtual ~ConstantNewStepper();

  virtual Real computeDT() override;

protected:
  const Real & _input_dt;

  const Real & _previous_stepper_dt;
};

#endif /* CONSTANTNEWSTEPPER_H */
