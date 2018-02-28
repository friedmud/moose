//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterVectorPostprocessors.h"

#include "Factory.h"

#include "CSVReader.h"
#include "ConstantVectorPostprocessor.h"
#include "Eigenvalues.h"
#include "ElementVariablesDifferenceMax.h"
#include "ElementsAlongLine.h"
#include "ElementsAlongPlane.h"
#include "IntersectionPointsAlongLine.h"
#include "LeastSquaresFit.h"
#include "LineFunctionSampler.h"
#include "LineMaterialRealSampler.h"
#include "LineValueSampler.h"
#include "MaterialVectorPostprocessor.h"
#include "NodalValueSampler.h"
#include "PointValueSampler.h"
#include "SideValueSampler.h"
#include "SphericalAverage.h"
#include "VectorOfPostprocessors.h"
#include "VolumeHistogram.h"

namespace Moose
{
/**
 * Called from Moose.C
 */
void
registerVectorPostprocessors(Factory & factory)
{
  registerVectorPostprocessor(CSVReader);
  registerVectorPostprocessor(ConstantVectorPostprocessor);
  registerVectorPostprocessor(Eigenvalues);
  registerVectorPostprocessor(ElementVariablesDifferenceMax);
  registerVectorPostprocessor(ElementsAlongLine);
  registerVectorPostprocessor(ElementsAlongPlane);
  registerVectorPostprocessor(IntersectionPointsAlongLine);
  registerVectorPostprocessor(LeastSquaresFit);
  registerVectorPostprocessor(LineFunctionSampler);
  registerVectorPostprocessor(LineMaterialRealSampler);
  registerVectorPostprocessor(LineValueSampler);
  registerVectorPostprocessor(MaterialVectorPostprocessor);
  registerVectorPostprocessor(NodalValueSampler);
  registerVectorPostprocessor(PointValueSampler);
  registerVectorPostprocessor(SideValueSampler);
  registerVectorPostprocessor(SphericalAverage);
  registerVectorPostprocessor(VectorOfPostprocessors);
  registerVectorPostprocessor(VolumeHistogram);
}
}
