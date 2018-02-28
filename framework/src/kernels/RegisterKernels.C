//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterKernels.h"

#include "Factory.h"

// Kernels
#include "ConservativeAdvection.h"
#include "TimeDerivative.h"
#include "CoupledTimeDerivative.h"
#include "MassLumpedTimeDerivative.h"
#include "Diffusion.h"
#include "AnisotropicDiffusion.h"
#include "CoupledForce.h"
#include "BodyForce.h"
#include "Reaction.h"
#include "MassEigenKernel.h"
#include "NullKernel.h"
#include "MaterialDerivativeTestKernel.h"
#include "MaterialDerivativeRankTwoTestKernel.h"
#include "MaterialDerivativeRankFourTestKernel.h"

// ScalarKernels
#include "ODETimeDerivative.h"
#include "CoupledODETimeDerivative.h"
#include "NodalEqualValueConstraint.h"
#include "ParsedODEKernel.h"
#include "QuotientScalarAux.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerKernels(Factory & factory)
{
  // Kernels
  registerKernel(TimeDerivative);
  registerKernel(ConservativeAdvection);
  registerKernel(CoupledTimeDerivative);
  registerKernel(MassLumpedTimeDerivative);
  registerKernel(Diffusion);
  registerKernel(AnisotropicDiffusion);
  registerKernel(CoupledForce);
  registerRenamedObject("UserForcingFunction", BodyForce, "04/01/2018 00:00");
  registerKernel(Reaction);
  registerKernel(MassEigenKernel);
  registerKernel(NullKernel);
  registerKernel(MaterialDerivativeTestKernel);
  registerKernel(MaterialDerivativeRankTwoTestKernel);
  registerKernel(MaterialDerivativeRankFourTestKernel);

  // Scalar Kernels
  registerScalarKernel(ODETimeDerivative);
  registerScalarKernel(CoupledODETimeDerivative);
  registerScalarKernel(NodalEqualValueConstraint);
  registerScalarKernel(ParsedODEKernel);
  registerScalarKernel(QuotientScalarAux);
}
}
