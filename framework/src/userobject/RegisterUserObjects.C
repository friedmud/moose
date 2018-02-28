//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterUserObjects.h"

#include "Factory.h"

#include "libmesh/libmesh_config.h"

#include "GeometrySphere.h"
#include "LayeredIntegral.h"
#include "LayeredAverage.h"
#include "LayeredSideIntegral.h"
#include "LayeredSideAverage.h"
#include "LayeredSideFluxAverage.h"
#include "NearestPointLayeredAverage.h"
#include "ElementIntegralVariableUserObject.h"
#include "NodalNormalsEvaluator.h"
#include "NodalNormalsCorner.h"
#include "NodalNormalsPreprocessor.h"
#include "SolutionUserObject.h"
#include "PerflogDumper.h"
#include "ElementQualityChecker.h"
#ifdef LIBMESH_HAVE_FPARSER
#include "Terminator.h"
#endif

namespace Moose
{
/**
 * Called from Moose.C
 */
void
registerUserObjects(Factory & factory)
{
  registerUserObject(GeometrySphere);
  registerUserObject(LayeredIntegral);
  registerUserObject(LayeredAverage);
  registerUserObject(LayeredSideIntegral);
  registerUserObject(LayeredSideAverage);
  registerUserObject(LayeredSideFluxAverage);
  registerUserObject(NearestPointLayeredAverage);
  registerUserObject(ElementIntegralVariableUserObject);
  registerUserObject(NodalNormalsPreprocessor);
  registerUserObject(NodalNormalsCorner);
  registerUserObject(NodalNormalsEvaluator);
  registerUserObject(SolutionUserObject);
  registerUserObject(PerflogDumper);
  registerUserObject(ElementQualityChecker);
#ifdef LIBMESH_HAVE_FPARSER
  registerUserObject(Terminator);
#endif
}
}
