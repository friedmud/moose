//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterTimeIntegrators.h"

#include "Factory.h"

#include "ImplicitEuler.h"
#include "BDF2.h"
#include "CrankNicolson.h"
#include "ExplicitEuler.h"
#include "ExplicitMidpoint.h"
#include "ExplicitTVDRK2.h"
#include "LStableDirk2.h"
#include "LStableDirk3.h"
#include "AStableDirk4.h"
#include "LStableDirk4.h"
#include "ImplicitMidpoint.h"
#include "Heun.h"
#include "Ralston.h"

namespace Moose
{
/**
 * Called from Moose.C
 */
void
registerTimeIntegrators(Factory & factory)
{
  registerTimeIntegrator(ImplicitEuler);
  registerTimeIntegrator(BDF2);
  registerTimeIntegrator(CrankNicolson);
  registerTimeIntegrator(ExplicitEuler);
  registerTimeIntegrator(ExplicitMidpoint);
  registerTimeIntegrator(ExplicitTVDRK2);
  registerTimeIntegrator(LStableDirk2);
  registerTimeIntegrator(LStableDirk3);
  registerTimeIntegrator(AStableDirk4);
  registerTimeIntegrator(LStableDirk4);
  registerTimeIntegrator(ImplicitMidpoint);
  registerTimeIntegrator(Heun);
  registerTimeIntegrator(Ralston);
}
}
