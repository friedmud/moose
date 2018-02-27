//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterDistributions.h"

#include "ConstantPointSource.h"
#include "FunctionDiracSource.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerDistributions(Factory & factory)
{
  registerDiracKernel(ConstantPointSource);
  registerDiracKernel(FunctionDiracSource);
}
}
