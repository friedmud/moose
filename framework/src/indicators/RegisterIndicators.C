//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterIndicators.h"

#include "Factory.h"

#include "AnalyticalIndicator.h"
#include "LaplacianJumpIndicator.h"
#include "GradientJumpIndicator.h"
#include "ValueJumpIndicator.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerIndicators(Factory & factory)
{
  registerIndicator(AnalyticalIndicator);
  registerIndicator(LaplacianJumpIndicator);
  registerIndicator(GradientJumpIndicator);
  registerIndicator(ValueJumpIndicator);
}
}
