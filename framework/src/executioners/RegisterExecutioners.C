//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterExecutioners.h"

#include "Factory.h"

#include "Steady.h"
#include "Transient.h"
#include "InversePowerMethod.h"
#include "NonlinearEigen.h"
#include "Eigenvalue.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerExecutioners(Factory & factory)
{
  registerExecutioner(Steady);
  registerExecutioner(Transient);
  registerExecutioner(InversePowerMethod);
  registerExecutioner(NonlinearEigen);
  registerExecutioner(Eigenvalue);
}
}
