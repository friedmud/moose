//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterNodalKernels.h"

#include "Factory.h"

#include "ConstantRate.h"
#include "TimeDerivativeNodalKernel.h"
#include "UserForcingFunctionNodalKernel.h"

namespace Moose
{
/**
 * Called from Moose.C
 */
void
registerNodalKernels(Factory & factory)
{
  registerNodalKernel(TimeDerivativeNodalKernel);
  registerNodalKernel(ConstantRate);
  registerNodalKernel(UserForcingFunctionNodalKernel);
}
}
