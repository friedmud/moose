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
class StepperInfo;

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

private:
  /// Information consumed by Steppers
  StepperInfo & _stepper_info;

protected:
  ////// Information from StepperInfo //////
  int & _step_count;
  Real & _time;
  std::deque<Real> & _dt;

  unsigned int & _nonlin_iters;
  unsigned int & _lin_iters;
  std::deque<bool> & _converged;
  std::deque<Real> & _solve_time_secs;

  std::unique_ptr<NumericVector<Number>> & _soln_nonlin;
  std::unique_ptr<NumericVector<Number>> & _soln_aux;
  std::unique_ptr<NumericVector<Number>> & _soln_predicted;

  bool & _snapshot;
  bool & _rewind;
  Real & _rewind_time;
};

#endif /* NEWSTEPPER_H */
