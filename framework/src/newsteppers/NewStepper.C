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
  params.addParam<bool>("reset_dt", false, "Use when restarting a calculation to force a change in dt.");

  params.registerBase("NewStepper");

  return params;
}

NewStepper::NewStepper(const InputParameters & parameters) :
    MooseObject(parameters),
    Restartable(parameters, "NewSteppers"),
    NewStepperInterface(this),
    _fe_problem(*parameters.getCheckedPointerParam<FEProblem *>("_fe_problem"))
{
}

NewStepper::~NewStepper()
{
}
