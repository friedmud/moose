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

#ifndef INITIALSTEPSNEWSTEPPER_H
#define INITIALSTEPSNEWSTEPPER_H

#include "NewStepper.h"

class InitialStepsNewStepper;

template<>
InputParameters validParams<InitialStepsNewStepper>();

/**
 * Cuts the incoming dt for the first N steps
 */
class InitialStepsNewStepper : public NewStepper
{
public:
  InitialStepsNewStepper(const InputParameters & parameters);
  virtual ~InitialStepsNewStepper();

  virtual Real computeDT() override;

  virtual Real computeFailedDT() override;

protected:
  const Real & _incoming_stepper_dt;

  const unsigned int & _n_steps;
};

#endif /* INITIALSTEPSNEWSTEPPER_H */
