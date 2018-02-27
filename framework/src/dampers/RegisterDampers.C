//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterDampers.h"

#include "ConstantDamper.h"
#include "MaxIncrement.h"
#include "BoundingValueNodalDamper.h"
#include "BoundingValueElementDamper.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerDampers(Factory & factory)
{
  registerDamper(ConstantDamper);
  registerDamper(MaxIncrement);
  registerDamper(BoundingValueNodalDamper);
  registerDamper(BoundingValueElementDamper);
}
}
