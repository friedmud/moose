//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterTransfers.h"

#include "Factory.h"

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_TRILINOS_HAVE_DTK
#include "MultiAppDTKUserObjectTransfer.h"
#include "MultiAppDTKInterpolationTransfer.h"
#endif
#include "MultiAppPostprocessorInterpolationTransfer.h"
#include "MultiAppVariableValueSampleTransfer.h"
#include "MultiAppVariableValueSamplePostprocessorTransfer.h"
#include "MultiAppMeshFunctionTransfer.h"
#include "MultiAppUserObjectTransfer.h"
#include "MultiAppNearestNodeTransfer.h"
#include "MultiAppCopyTransfer.h"
#include "MultiAppInterpolationTransfer.h"
#include "MultiAppPostprocessorTransfer.h"
#include "MultiAppProjectionTransfer.h"
#include "MultiAppPostprocessorToAuxScalarTransfer.h"
#include "MultiAppScalarToAuxScalarTransfer.h"
#include "MultiAppVectorPostprocessorTransfer.h"

namespace Moose
{
/**
 * Called from Moose.C
 */
void
registerTransfers(Factory & factory)
{
#ifdef LIBMESH_TRILINOS_HAVE_DTK
  registerTransfer(MultiAppDTKUserObjectTransfer);
  registerTransfer(MultiAppDTKInterpolationTransfer);
#endif
  registerTransfer(MultiAppPostprocessorInterpolationTransfer);
  registerTransfer(MultiAppVariableValueSampleTransfer);
  registerTransfer(MultiAppVariableValueSamplePostprocessorTransfer);
  registerTransfer(MultiAppMeshFunctionTransfer);
  registerTransfer(MultiAppUserObjectTransfer);
  registerTransfer(MultiAppNearestNodeTransfer);
  registerTransfer(MultiAppCopyTransfer);
  registerTransfer(MultiAppInterpolationTransfer);
  registerTransfer(MultiAppPostprocessorTransfer);
  registerTransfer(MultiAppProjectionTransfer);
  registerTransfer(MultiAppPostprocessorToAuxScalarTransfer);
  registerTransfer(MultiAppScalarToAuxScalarTransfer);
  registerTransfer(MultiAppVectorPostprocessorTransfer);
}
}
