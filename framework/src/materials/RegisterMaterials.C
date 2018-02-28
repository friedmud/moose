//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterMaterials.h"

#include "Factory.h"

#include "DerivativeParsedMaterial.h"
#include "DerivativeSumMaterial.h"
#include "GenericConstantMaterial.h"
#include "GenericConstantRankTwoTensor.h"
#include "GenericFunctionMaterial.h"
#include "ParsedMaterial.h"
#include "PiecewiseLinearInterpolationMaterial.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerMaterials(Factory & factory)
{
  registerMaterial(DerivativeParsedMaterial);
  registerMaterial(DerivativeSumMaterial);
  registerMaterial(GenericConstantMaterial);
  registerMaterial(GenericConstantRankTwoTensor);
  registerMaterial(GenericFunctionMaterial);
  registerMaterial(ParsedMaterial);
  registerMaterial(PiecewiseLinearInterpolationMaterial);
}
}
