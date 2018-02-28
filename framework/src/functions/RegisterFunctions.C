//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterFunctions.h"

#include "Factory.h"

#include "Axisymmetric2D3DSolutionFunction.h"
#include "ConstantFunction.h"
#include "CompositeFunction.h"
#include "MooseParsedFunction.h"
#include "MooseParsedVectorFunction.h"
#include "MooseParsedGradFunction.h"
#include "PiecewiseConstant.h"
#include "PiecewiseLinear.h"
#include "SolutionFunction.h"
#include "PiecewiseBilinear.h"
#include "SplineFunction.h"
#include "BicubicSplineFunction.h"
#include "PiecewiseMultilinear.h"
#include "LinearCombinationFunction.h"
#include "ImageFunction.h"
#include "VectorPostprocessorFunction.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerFunctions(Factory & factory)
{
  registerFunction(Axisymmetric2D3DSolutionFunction);
  registerFunction(ConstantFunction);
  registerFunction(CompositeFunction);
  registerNamedFunction(MooseParsedFunction, "ParsedFunction");
  registerNamedFunction(MooseParsedGradFunction, "ParsedGradFunction");
  registerNamedFunction(MooseParsedVectorFunction, "ParsedVectorFunction");
  registerFunction(PiecewiseConstant);
  registerFunction(PiecewiseLinear);
  registerFunction(SolutionFunction);
  registerFunction(PiecewiseBilinear);
  registerFunction(SplineFunction);
  registerFunction(BicubicSplineFunction);
  registerFunction(PiecewiseMultilinear);
  registerFunction(LinearCombinationFunction);
  registerFunction(ImageFunction);
  registerFunction(VectorPostprocessorFunction);
}
}
