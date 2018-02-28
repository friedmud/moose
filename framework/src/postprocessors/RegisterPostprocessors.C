//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterPostprocessors.h"

#include "Factory.h"

#include "AverageElementSize.h"
#include "AverageNodalVariableValue.h"
#include "CumulativeValuePostprocessor.h"
#include "ChangeOverTimePostprocessor.h"
#include "ChangeOverTimestepPostprocessor.h"
#include "NodalSum.h"
#include "ElementAverageValue.h"
#include "ElementAverageTimeDerivative.h"
#include "ElementW1pError.h"
#include "ElementH1Error.h"
#include "ElementH1SemiError.h"
#include "ElementIntegralVariablePostprocessor.h"
#include "ElementIntegralMaterialProperty.h"
#include "ElementL2Error.h"
#include "ElementVectorL2Error.h"
#include "EmptyPostprocessor.h"
#include "FindValueOnLine.h"
#include "FunctionValuePostprocessor.h"
#include "NodalVariableValue.h"
#include "NumDOFs.h"
#include "TimestepSize.h"
#include "PerformanceData.h"
#include "MemoryUsage.h"
#include "NumElems.h"
#include "NumNodes.h"
#include "NumNonlinearIterations.h"
#include "NumLinearIterations.h"
#include "Residual.h"
#include "ScalarVariable.h"
#include "NumVars.h"
#include "NumResidualEvaluations.h"
#include "Receiver.h"
#include "SideAverageValue.h"
#include "SideFluxIntegral.h"
#include "SideFluxAverage.h"
#include "SideIntegralVariablePostprocessor.h"
#include "NodalMaxValue.h"
#include "NodalProxyMaxValue.h"
#include "ScalarL2Error.h"
#include "ElementalVariableValue.h"
#include "ElementL2Norm.h"
#include "NodalL2Norm.h"
#include "NodalL2Error.h"
#include "TotalVariableValue.h"
#include "VolumePostprocessor.h"
#include "AreaPostprocessor.h"
#include "PointValue.h"
#include "NodalExtremeValue.h"
#include "ElementExtremeValue.h"
#include "DifferencePostprocessor.h"
#include "RelativeDifferencePostprocessor.h"
#include "ScalePostprocessor.h"
#include "LinearCombinationPostprocessor.h"
#include "NumPicardIterations.h"
#include "FunctionSideIntegral.h"
#include "ExecutionerAttributeReporter.h"
#include "PercentChangePostprocessor.h"
#include "ElementL2Difference.h"
#include "TimeExtremeValue.h"
#include "RelativeSolutionDifferenceNorm.h"
#include "AxisymmetricCenterlineAverageValue.h"
#include "VariableInnerProduct.h"
#include "VariableResidual.h"

namespace Moose
{
/**
 * Called from Moose.C
 */
void
registerPostprocessors(Factory & factory)
{
  registerPostprocessor(AverageElementSize);
  registerPostprocessor(AverageNodalVariableValue);
  registerPostprocessor(CumulativeValuePostprocessor);
  registerPostprocessor(ChangeOverTimePostprocessor);
  registerPostprocessor(ChangeOverTimestepPostprocessor);
  registerPostprocessor(NodalSum);
  registerPostprocessor(ElementAverageValue);
  registerPostprocessor(ElementAverageTimeDerivative);
  registerPostprocessor(ElementW1pError);
  registerPostprocessor(ElementH1Error);
  registerPostprocessor(ElementH1SemiError);
  registerPostprocessor(ElementIntegralVariablePostprocessor);
  registerPostprocessor(ElementIntegralMaterialProperty);
  registerPostprocessor(ElementL2Error);
  registerPostprocessor(ElementVectorL2Error);
  registerPostprocessor(ScalarL2Error);
  registerPostprocessor(EmptyPostprocessor);
  registerPostprocessor(FindValueOnLine);
  registerPostprocessor(NodalVariableValue);
  registerPostprocessor(NumDOFs);
  registerPostprocessor(TimestepSize);
  registerPostprocessor(PerformanceData);
  registerPostprocessor(MemoryUsage);
  registerPostprocessor(NumElems);
  registerPostprocessor(NumNodes);
  registerPostprocessor(NumNonlinearIterations);
  registerPostprocessor(NumLinearIterations);
  registerPostprocessor(Residual);
  registerPostprocessor(ScalarVariable);
  registerPostprocessor(NumVars);
  registerPostprocessor(NumResidualEvaluations);
  registerPostprocessor(Receiver);
  registerPostprocessor(SideAverageValue);
  registerPostprocessor(SideFluxIntegral);
  registerPostprocessor(SideFluxAverage);
  registerPostprocessor(SideIntegralVariablePostprocessor);
  registerPostprocessor(NodalMaxValue);
  registerPostprocessor(NodalProxyMaxValue);
  registerPostprocessor(ElementalVariableValue);
  registerPostprocessor(ElementL2Norm);
  registerPostprocessor(NodalL2Norm);
  registerPostprocessor(NodalL2Error);
  registerPostprocessor(TotalVariableValue);
  registerPostprocessor(VolumePostprocessor);
  registerPostprocessor(AreaPostprocessor);
  registerPostprocessor(PointValue);
  registerPostprocessor(NodalExtremeValue);
  registerPostprocessor(ElementExtremeValue);
  registerPostprocessor(DifferencePostprocessor);
  registerPostprocessor(RelativeDifferencePostprocessor);
  registerPostprocessor(ScalePostprocessor);
  registerPostprocessor(LinearCombinationPostprocessor);
  registerPostprocessor(FunctionValuePostprocessor);
  registerPostprocessor(NumPicardIterations);
  registerPostprocessor(FunctionSideIntegral);
  registerPostprocessor(ExecutionerAttributeReporter);
  registerPostprocessor(PercentChangePostprocessor);
  registerPostprocessor(ElementL2Difference);
  registerPostprocessor(TimeExtremeValue);
  registerPostprocessor(RelativeSolutionDifferenceNorm);
  registerPostprocessor(AxisymmetricCenterlineAverageValue);
  registerPostprocessor(VariableInnerProduct);
  registerPostprocessor(VariableResidual);
}
}
