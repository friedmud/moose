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

#include "CuttingGrowingNewStepper.h"
#include "FEProblem.h"
#include "Transient.h"
#include "MooseApp.h"

template<>
InputParameters validParams<CuttingGrowingNewStepper>()
{
  InputParameters params = validParams<NewStepper>();

  params.addRequiredParam<Real>("dt", "The requested dt.  The dt will be cut when a solve fails and will attempt to regrow back to this value when solves converge.");

  params.addRangeCheckedParam<Real>("growth_factor", 2, "growth_factor>=1", "Maximum ratio of new to previous timestep sizes following a step that required the time step to be cut due to a failed solve.");

  return params;
}

CuttingGrowingNewStepper::CuttingGrowingNewStepper(const InputParameters & parameters) :
    NewStepper(parameters),
    _input_dt(getParam<Real>("dt"))
{
  // Start with the old dt
  {
    auto params = _factory.getValidParams("OldDTNewStepper");
    _fe_problem->addNewStepper("OldDTNewStepper", name() + "_old_dt", params);
  }

  // Now multiply it by the growth factor
  {
    auto params = _factory.getValidParams("MultiplierNewStepper");
    params.set<Real>("multiplier") = _growth_factor;
    params.set<StepperName>("previous_stepper") = name() + "_old_dt";
    _fe_problem->addNewStepper("MultiplierNewStepper", name() + "_growth_mult", params);
  }

  // Also multiply it by 0.5
  {
    auto params = _factory.getValidParams("MultiplierNewStepper");
    params.set<Real>("multiplier") = 0.5;
    params.set<StepperName>("previous_stepper") = name() + "_old_dt";
    _fe_problem->addNewStepper("MultiplierNewStepper", name() + "_cut_mult", params);
  }

  // Choose one of those based on whether or not it converged
  {
    auto params = _factory.getValidParams("ConvergedNewStepper");
    params.set<StepperName>("converged") = name() + "_growth_mult";
    params.set<StepperName>("unconverged") = name() + "_cut_mult";
    _fe_problem->addNewStepper("ConvergedNewStepper", name() + "_converged1", params);
  }

  // Create a ConstantNewStepper with the dt passed into this object
  {
    auto params = _factory.getValidParams("ConstantNewStepper");
    params.set<Real>("dt") = _input_dt;
    _fe_problem->addNewStepper("ConstantNewStepper", name() + "_constant", params);
  }

  // Choose the minimum from the converged above and the specified dt
  {
    auto params = _factory.getValidParams("MinumumNewStepper");
    params.set<StepperName>("previous_stepper1") = name() + "_converged1";
    params.set<StepperName>("previous_stepper2") = name() + "_constant";
    _fe_problem->addNewStepper("MinimumNewStepper", name() + "_min1", params);
  }

  // Honor startup steps
  {
    auto params = _factory.getValidParams("StartupNewStepper");
    params.set<StepperName>("previous_stepper") = name() + "_min1";
    _fe_problem->addNewStepper("StartupNewStepper", name() + "_startup", params);
  }

  // Finally, multiply... I don't like this...
  {
  }
}

Real
CuttingGrowingNewStepper::computeDT()
{
  return _input_dt;
}

CuttingGrowingNewStepper::~CuttingGrowingNewStepper()
{
}
