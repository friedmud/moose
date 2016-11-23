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

#include "ConstantNewStepper.h"
#include "FEProblem.h"
#include "Transient.h"
#include "MooseApp.h"

template<>
InputParameters validParams<ConstantNewStepper>()
{
  InputParameters params = validParams<NewStepper>();

  params.addRequiredParam<Real>("dt", "The dt to maintain");
  params.addParam<StepperName>("previous_stepper", "", "The name of the previous Stepper to pull dt from");

  return params;
}

ConstantNewStepper::ConstantNewStepper(const InputParameters & parameters) :
    NewStepper(parameters),
    _input_dt(getParam<Real>("dt")),
    _previous_stepper_dt(getStepperDT("previous_stepper"))
{
}

Real
ConstantNewStepper::computeDT()
{
  return _input_dt;
}

ConstantNewStepper::~ConstantNewStepper()
{
}
