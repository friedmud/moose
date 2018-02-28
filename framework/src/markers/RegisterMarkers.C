//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterMarkers.h"

#include "Factory.h"

#include "ErrorToleranceMarker.h"
#include "ErrorFractionMarker.h"
#include "UniformMarker.h"
#include "BoxMarker.h"
#include "ComboMarker.h"
#include "ValueThresholdMarker.h"
#include "ValueRangeMarker.h"
#include "OrientedBoxMarker.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerMarkers(Factory & factory)
{
  registerMarker(ErrorToleranceMarker);
  registerMarker(ErrorFractionMarker);
  registerMarker(UniformMarker);
  registerMarker(BoxMarker);
  registerMarker(OrientedBoxMarker);
  registerMarker(ComboMarker);
  registerMarker(ValueThresholdMarker);
  registerMarker(ValueRangeMarker);
}
}
