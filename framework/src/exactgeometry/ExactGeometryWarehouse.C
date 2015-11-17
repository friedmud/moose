/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ExactGeometryWarehouse.h"
#include "ExactGeometry.h"


ExactGeometryWarehouse::ExactGeometryWarehouse()
{
}

ExactGeometryWarehouse::~ExactGeometryWarehouse()
{
}

void
ExactGeometryWarehouse::initialSetup()
{
  for (std::vector<ExactGeometry *>::const_iterator i = all().begin(); i != all().end(); ++i)
    (*i)->initialSetup();
}

void
ExactGeometryWarehouse::timestepSetup()
{
  for (std::vector<ExactGeometry *>::const_iterator i = all().begin(); i != all().end(); ++i)
    (*i)->timestepSetup();
}

void
ExactGeometryWarehouse::residualSetup()
{
  for (std::vector<ExactGeometry *>::const_iterator i = all().begin(); i != all().end(); ++i)
    (*i)->residualSetup();
}

void
ExactGeometryWarehouse::jacobianSetup()
{
  for (std::vector<ExactGeometry *>::const_iterator i = all().begin(); i != all().end(); ++i)
    (*i)->jacobianSetup();
}

void
ExactGeometryWarehouse::addExactGeometry(MooseSharedPointer<ExactGeometry> exact_geometry)
{
  // Make certain that the ExactGeometry is valid
  mooseAssert(exact_geometry, "ExactGeometry is NULL");

  _all_objects.push_back(exact_geometry.get());

  // Add to elemental/nodal storage
  _all_exact_geometries.push_back(exact_geometry);

  const std::set<BoundaryID> & boundaries = exact_geometry->boundaryIDs();

  for (std::set<BoundaryID>::const_iterator bit = boundaries.begin(); bit != boundaries.end(); ++bit)
  {
    BoundaryID bnd_id = *bit;

    _active_boundary_exact_geometries[bnd_id].push_back(exact_geometry);
  }
}
