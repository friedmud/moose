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

#ifndef NEWSTEPPER_H
#define NEWSTEPPER_H

#include "MooseObject.h"
#include "Restartable.h"
#include "NewStepperInterface.h"

class NewStepper;
class FEProblem;
class Transient;
class StepperBlock;

template<>
InputParameters validParams<NewStepper>();

/**
 * Base class for time stepping
 */
class NewStepper :
  public MooseObject,
  public Restartable,
  public NewStepperInterface
{
public:
  NewStepper(const InputParameters & parameters);
  virtual ~NewStepper();

  /**
   * Compute the size of the next timestep
   */
  virtual Real computeDT() = 0;

protected:
  FEProblem & _fe_problem;
};

#endif /* NEWSTEPPER_H */
