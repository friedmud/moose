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

#include "NewStepper.h"
#include "FEProblem.h"
#include "Transient.h"
#include "MooseApp.h"

template<>
InputParameters validParams<NewStepper>()
{
  InputParameters params = validParams<MooseObject>();

  params.registerBase("NewStepper");

  return params;
}

NewStepper::NewStepper(const InputParameters & parameters) :
    MooseObject(parameters),
    Restartable(parameters, "NewSteppers"),
    NewStepperInterface(this),
    _fe_problem(*parameters.getCheckedPointerParam<FEProblem *>("_fe_problem")),
    _factory(_app.getFactory()),
    _stepper_info(_fe_problem.getStepperInfo()),
    _step_count(_stepper_info._step_count),
    _time(_stepper_info._time),
    _dt(_stepper_info._dt),
    _nonlin_iters(_stepper_info._nonlin_iters),
    _lin_iters(_stepper_info._lin_iters),
    _converged(_stepper_info._converged),
    _solve_time_secs(_stepper_info._solve_time_secs),
    _soln_nonlin(_stepper_info._soln_nonlin),
    _soln_aux(_stepper_info._soln_aux),
    _soln_predicted(_stepper_info._soln_predicted),
    _snapshot(_stepper_info._snapshot),
    _rewind(_stepper_info._rewind),
    _rewind_time(_stepper_info._rewind_time)
{
}

NewStepper::~NewStepper()
{
}
