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

#include "NewStepperInterface.h"
#include "FEProblem.h"
#include "MooseTypes.h"
#include "MooseObject.h"

NewStepperInterface::NewStepperInterface(const MooseObject * moose_object) :
    _si_params(moose_object->parameters()),
    _si_feproblem(*_si_params.getCheckedPointerParam<FEProblem *>("_fe_problem",
                                                                  "Missing FEProblem Pointer in NewStepperInterface!"))
{
}

const Real &
NewStepperInterface::getStepperDT(const std::string & name)
{
  if (!hasStepper(name))
    mooseError("Unknown Stepper parameter name: " << name);

  return _si_feproblem.getStepperDT(_si_params.get<StepperName>(name));
}

const Real &
NewStepperInterface::getStepperDTByName(const StepperName & name)
{
  if (!hasStepperByName(name))
    mooseError("Unknown Stepper: " << name);

  return _si_feproblem.getStepperDT(name);
}

bool
NewStepperInterface::hasStepper(const std::string & name) const
{
  return _si_feproblem.hasPostprocessor(_si_params.get<PostprocessorName>(name));
}

bool
NewStepperInterface::hasStepperByName(const PostprocessorName & name) const
{
  return _si_feproblem.hasPostprocessor(name);
}
