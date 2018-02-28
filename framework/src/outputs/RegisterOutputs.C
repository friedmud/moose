//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterOutputs.h"

#include "Factory.h"

#include "libmesh/petsc_macro.h"
#include "libmesh/libmesh_config.h"

#include "Moose.h"

#ifdef LIBMESH_HAVE_EXODUS_API
#include "Exodus.h"
#endif
#include "Nemesis.h"
#include "Console.h"
#include "CSV.h"
#include "VTKOutput.h"
#include "Checkpoint.h"
#include "XDA.h"
#include "GMVOutput.h"
#include "Tecplot.h"
#include "Gnuplot.h"
#include "SolutionHistory.h"
#include "MaterialPropertyDebugOutput.h"
#include "VariableResidualNormsDebugOutput.h"
#include "TopResidualDebugOutput.h"
#include "DOFMapOutput.h"
#include "ControlOutput.h"

namespace Moose
{
/**
 * Called from Moose.C to register Actions
 */
void
registerOutputs(Factory & factory)
{
#ifdef LIBMESH_HAVE_EXODUS_API
  registerOutput(Exodus);
#endif
#ifdef LIBMESH_HAVE_NEMESIS_API
  registerOutput(Nemesis);
#endif
  registerOutput(Console);
  registerOutput(CSV);
#ifdef LIBMESH_HAVE_VTK
  registerNamedOutput(VTKOutput, "VTK");
#endif
  registerOutput(Checkpoint);
  registerNamedOutput(XDA, "XDR");
  registerOutput(XDA);
  registerNamedOutput(GMVOutput, "GMV");
  registerOutput(Tecplot);
  registerOutput(Gnuplot);
  registerOutput(SolutionHistory);
  registerOutput(MaterialPropertyDebugOutput);
  registerOutput(VariableResidualNormsDebugOutput);
  registerOutput(TopResidualDebugOutput);
  registerNamedOutput(DOFMapOutput, "DOFMap");
  registerOutput(ControlOutput);
}
}
