//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterPreconditioners.h"

#include "Factory.h"

#include "PhysicsBasedPreconditioner.h"
#include "FiniteDifferencePreconditioner.h"
#include "SingleMatrixPreconditioner.h"
#include "FieldSplitPreconditioner.h"
#include "Split.h"

namespace Moose
{
/**
 * Called from Moose.C
 */
void
registerPreconditioners(Factory & factory)
{
  registerNamedPreconditioner(PhysicsBasedPreconditioner, "PBP");
  registerNamedPreconditioner(FiniteDifferencePreconditioner, "FDP");
  registerNamedPreconditioner(SingleMatrixPreconditioner, "SMP");
#if defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3, 3, 0)
  registerNamedPreconditioner(FieldSplitPreconditioner, "FSP");
#endif
}
}
