//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterMeshModifiers.h"

#include "Factory.h"

#include "MeshExtruder.h"
#include "SideSetsFromPoints.h"
#include "SideSetsFromNormals.h"
#include "AddExtraNodeset.h"
#include "BoundingBoxNodeSet.h"
#include "Transform.h"
#include "SideSetsAroundSubdomain.h"
#include "SideSetsBetweenSubdomains.h"
#include "AddAllSideSetsByNormals.h"
#include "SubdomainBoundingBox.h"
#include "OrientedSubdomainBoundingBox.h"
#include "RenameBlock.h"
#include "AssignElementSubdomainID.h"
#include "ImageSubdomain.h"
#include "BlockDeleter.h"
#include "ParsedSubdomainMeshModifier.h"
#include "BreakBoundaryOnSubdomain.h"
#include "ParsedAddSideset.h"
#include "AssignSubdomainID.h"
#include "MeshSideSet.h"
#include "AddSideSetsFromBoundingBox.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerMeshModifiers(Factory & factory)
{
  registerMeshModifier(MeshExtruder);
  registerMeshModifier(SideSetsFromPoints);
  registerMeshModifier(SideSetsFromNormals);
  registerMeshModifier(AddExtraNodeset);
  registerMeshModifier(BoundingBoxNodeSet);
  registerMeshModifier(Transform);
  registerMeshModifier(SideSetsAroundSubdomain);
  registerMeshModifier(SideSetsBetweenSubdomains);
  registerMeshModifier(AddAllSideSetsByNormals);
  registerMeshModifier(SubdomainBoundingBox);
  registerMeshModifier(OrientedSubdomainBoundingBox);
  registerMeshModifier(RenameBlock);
  registerMeshModifier(AssignElementSubdomainID);
  registerMeshModifier(ImageSubdomain);
  registerMeshModifier(BlockDeleter);
  registerMeshModifier(ParsedSubdomainMeshModifier);
  registerMeshModifier(BreakBoundaryOnSubdomain);
  registerMeshModifier(ParsedAddSideset);
  registerMeshModifier(AssignSubdomainID);
  registerMeshModifier(MeshSideSet);
  registerMeshModifier(AddSideSetsFromBoundingBox);
}
}
