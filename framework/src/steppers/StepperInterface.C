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

#include "StepperInterface.h"
#include "FEProblem.h"
#include "MooseTypes.h"
#include "MooseObject.h"

StepperInterface::StepperInterface(const MooseObject * moose_object) :
    DependencyResolverInterface(),
    _si_name({moose_object->parameters().get<StepperName>("_output_name") != "" ? moose_object->parameters().get<StepperName>("_output_name") : moose_object->parameters().get<std::string>("_object_name")}),
    _si_params(moose_object->parameters()),
    _si_feproblem(*_si_params.getCheckedPointerParam<FEProblem *>("_fe_problem",
                                                                  "Missing FEProblem Pointer in StepperInterface!"))
{
}

const Real &
StepperInterface::getStepperDT(const std::string & name)
{
  std::cout<<"StepperInterface::getStepperDT: "<<name<<std::endl;
  std::cout<<"StepperInterface::getStepperDT: "<<_si_params.get<StepperName>(name)<<std::endl;

  _depend_steppers.insert(_si_params.get<StepperName>(name));

  return _si_feproblem.getStepperDT(_si_params.get<StepperName>(name));
}

const Real &
StepperInterface::getStepperDTByName(const StepperName & name)
{
  _depend_steppers.insert(name);

  return _si_feproblem.getStepperDT(name);
}

const std::set<std::string> &
StepperInterface::getRequestedItems()
{
  return _depend_steppers;
}

const std::set<std::string> &
StepperInterface::getSuppliedItems()
{
  return _si_name;
}

void
StepperInterface::setSuppliedItemName(const StepperName & item_name)
{
  _si_name = {static_cast<std::string>(item_name)};
}

bool
StepperInterface::hasStepper(const std::string & name) const
{
  return _si_feproblem.hasPostprocessor(_si_params.get<PostprocessorName>(name));
}

bool
StepperInterface::hasStepperByName(const PostprocessorName & name) const
{
  return _si_feproblem.hasPostprocessor(name);
}
