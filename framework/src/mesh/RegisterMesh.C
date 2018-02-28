//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterMesh.h"

#include "Factory.h"

#include "FileMesh.h"
#include "GeneratedMesh.h"
#include "TiledMesh.h"
#include "ImageMesh.h"
#include "PatternedMesh.h"
#include "StitchedMesh.h"
#include "AnnularMesh.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerMeshObjects(Factory & factory)
{
  registerMesh(FileMesh);
  registerMesh(GeneratedMesh);
  registerMesh(TiledMesh);
  registerMesh(ImageMesh);
  registerMesh(PatternedMesh);
  registerMesh(StitchedMesh);
  registerMesh(AnnularMesh);
}
}
