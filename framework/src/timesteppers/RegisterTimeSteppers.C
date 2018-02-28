//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterTimeSteppers.h"

#include "Factory.h"

#include "ConstantDT.h"
#include "LogConstantDT.h"
#include "FunctionDT.h"
#include "TimeSequenceStepper.h"
#include "ExodusTimeSequenceStepper.h"
#include "CSVTimeSequenceStepper.h"
#include "IterationAdaptiveDT.h"
#include "SolutionTimeAdaptiveDT.h"
#include "DT2.h"
#include "PostprocessorDT.h"
#include "AB2PredictorCorrector.h"

namespace Moose
{
/**
 * Called from Moose.C
 */
void
registerTimeSteppers(Factory & factory)
{
  registerTimeStepper(ConstantDT);
  registerTimeStepper(LogConstantDT);
  registerTimeStepper(FunctionDT);
  registerTimeStepper(TimeSequenceStepper);
  registerTimeStepper(ExodusTimeSequenceStepper);
  registerTimeStepper(CSVTimeSequenceStepper);
  registerTimeStepper(IterationAdaptiveDT);
  registerTimeStepper(SolutionTimeAdaptiveDT);
  registerTimeStepper(DT2);
  registerTimeStepper(PostprocessorDT);
  registerTimeStepper(AB2PredictorCorrector);
}
}
