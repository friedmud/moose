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

#ifndef BETTERCUTTINGGROWINGNEWSTEPPER_H
#define BETTERCUTTINGGROWINGNEWSTEPPER_H

#include "NewStepper.h"

class BetterCuttingGrowingNewStepper;

template<>
InputParameters validParams<BetterCuttingGrowingNewStepper>();

/**
 * Essentially takes the place of what used to be ConstantDT
 * The name was so terribly bad for that old one that we gotta get rid of it.
 * We might cross-register this object with ConstantDT once it's deployed though
 * and register it as a deprecated name.
 *
 * This one cuts the timestep by half when a solve fails and regrows by
 * a growth_factor when there is convergence.
 */
class BetterCuttingGrowingNewStepper : public NewStepper
{
public:
  BetterCuttingGrowingNewStepper(const InputParameters & parameters);
  virtual ~BetterCuttingGrowingNewStepper();

  virtual Real computeDT() override;

  virtual Real computeFailedDT() override;

protected:
  const Real & _input_dt;

  /// The factor to grow dt by
  const Real & _growth_factor;

  /// The previous dt _this_ particular Stepper computed
  Real & _prev_dt;
};

#endif /* BETTERCUTTINGGROWINGNEWSTEPPER_H */
