//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterICs.h"

#include "Factory.h"

#include "ConstantIC.h"
#include "BoundingBoxIC.h"
#include "FunctionIC.h"
#include "RandomIC.h"
#include "ScalarConstantIC.h"
#include "ScalarComponentIC.h"
#include "FunctionScalarIC.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerICs(Factory & factory)
{
  registerInitialCondition(ConstantIC);
  registerInitialCondition(BoundingBoxIC);
  registerInitialCondition(FunctionIC);
  registerInitialCondition(RandomIC);
  registerInitialCondition(ScalarConstantIC);
  registerInitialCondition(ScalarComponentIC);
  registerInitialCondition(FunctionScalarIC);
}
}
