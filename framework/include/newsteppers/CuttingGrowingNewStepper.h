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

#ifndef CUTTINGGROWINGNEWSTEPPER_H
#define CUTTINGGROWINGNEWSTEPPER_H

#include "NewStepper.h"

class CuttingGrowingNewStepper;

template<>
InputParameters validParams<CuttingGrowingNewStepper>();

/**
 * Essentially takes the place of what used to be ConstantDT
 * The name was so terribly bad for that old one that we gotta get rid of it.
 * We might cross-register this object with ConstantDT once it's deployed though
 * and register it as a deprecated name.
 *
 * This one cuts the timestep by half when a solve fails and regrows by
 * a growth_factor when there is convergence.
 */
class CuttingGrowingNewStepper : public NewStepper
{
public:
  CuttingGrowingNewStepper(const InputParameters & parameters);
  virtual ~CuttingGrowingNewStepper();

  virtual Real computeDT() override;

protected:
  const Real & _input_dt;

  const Real & _previous_stepper_dt;
};

#endif /* CUTTINGGROWINGNEWSTEPPER_H */
