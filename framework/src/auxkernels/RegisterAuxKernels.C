//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterAuxKernels.h"

#include "Factory.h"

#include "ConstantAux.h"
#include "FunctionAux.h"
#include "NearestNodeDistanceAux.h"
#include "NearestNodeValueAux.h"
#include "PenetrationAux.h"
#include "ProcessorIDAux.h"
#include "SelfAux.h"
#include "GapValueAux.h"
#include "MaterialRealAux.h"
#include "MaterialRealVectorValueAux.h"
#include "MaterialRealTensorValueAux.h"
#include "MaterialStdVectorAux.h"
#include "MaterialRealDenseMatrixAux.h"
#include "MaterialStdVectorRealGradientAux.h"
#include "DebugResidualAux.h"
#include "BoundsAux.h"
#include "SpatialUserObjectAux.h"
#include "SolutionAux.h"
#include "VectorMagnitudeAux.h"
#include "ConstantScalarAux.h"
#include "QuotientAux.h"
#include "NormalizationAux.h"
#include "VariableGradientComponent.h"
#include "ParsedAux.h"
#include "VariableTimeIntegrationAux.h"
#include "ElementLengthAux.h"
#include "ElementLpNormAux.h"
#include "ElementL2ErrorFunctionAux.h"
#include "ElementH1ErrorFunctionAux.h"
#include "DiffusionFluxAux.h"
#include "FunctionScalarAux.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerAuxKernels(Factory & factory)
{
  registerAux(ConstantAux);
  registerAux(FunctionAux);
  registerAux(NearestNodeDistanceAux);
  registerAux(NearestNodeValueAux);
  registerAux(PenetrationAux);
  registerAux(ProcessorIDAux);
  registerAux(SelfAux);
  registerAux(GapValueAux);
  registerAux(MaterialRealAux);
  registerAux(MaterialRealVectorValueAux);
  registerAux(MaterialRealTensorValueAux);
  registerAux(MaterialStdVectorAux);
  registerAux(MaterialRealDenseMatrixAux);
  registerAux(MaterialStdVectorRealGradientAux);
  registerAux(DebugResidualAux);
  registerAux(BoundsAux);
  registerAux(SpatialUserObjectAux);
  registerAux(SolutionAux);
  registerAux(VectorMagnitudeAux);
  registerAux(ConstantScalarAux);
  registerAux(QuotientAux);
  registerAux(NormalizationAux);
  registerAux(FunctionScalarAux);
  registerAux(VariableGradientComponent);
  registerAux(ParsedAux);
  registerAux(VariableTimeIntegrationAux);
  registerAux(ElementLengthAux);
  registerAux(ElementLpNormAux);
  registerAux(ElementL2ErrorFunctionAux);
  registerAux(ElementH1ErrorFunctionAux);
  registerAux(DiffusionFluxAux);
}
}
