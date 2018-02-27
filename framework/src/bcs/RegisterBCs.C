//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterBCs.h"

#include "ConvectiveFluxBC.h"
#include "DirichletBC.h"
#include "PenaltyDirichletBC.h"
#include "PresetBC.h"
#include "NeumannBC.h"
#include "PostprocessorNeumannBC.h"
#include "FunctionDirichletBC.h"
#include "FunctionPenaltyDirichletBC.h"
#include "FunctionPresetBC.h"
#include "FunctionNeumannBC.h"
#include "MatchedValueBC.h"
#include "VacuumBC.h"
#include "SinDirichletBC.h"
#include "SinNeumannBC.h"
#include "VectorNeumannBC.h"
#include "WeakGradientBC.h"
#include "DiffusionFluxBC.h"
#include "PostprocessorDirichletBC.h"
#include "OneDEqualValueConstraintBC.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerBCs(Factory & factory)
{
  registerBoundaryCondition(ConvectiveFluxBC);
  registerBoundaryCondition(DirichletBC);
  registerBoundaryCondition(PenaltyDirichletBC);
  registerBoundaryCondition(PresetBC);
  registerBoundaryCondition(NeumannBC);
  registerBoundaryCondition(PostprocessorNeumannBC);
  registerBoundaryCondition(FunctionDirichletBC);
  registerBoundaryCondition(FunctionPenaltyDirichletBC);
  registerBoundaryCondition(FunctionPresetBC);
  registerBoundaryCondition(FunctionNeumannBC);
  registerBoundaryCondition(MatchedValueBC);
  registerBoundaryCondition(VacuumBC);
  registerBoundaryCondition(SinDirichletBC);
  registerBoundaryCondition(SinNeumannBC);
  registerBoundaryCondition(VectorNeumannBC);
  registerBoundaryCondition(WeakGradientBC);
  registerBoundaryCondition(DiffusionFluxBC);
  registerBoundaryCondition(PostprocessorDirichletBC);
  registerBoundaryCondition(OneDEqualValueConstraintBC);
}
}
