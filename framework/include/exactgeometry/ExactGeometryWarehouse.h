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

#ifndef EXACTGEOMETRYWAREHOUSE_H
#define EXACTGEOMETRYWAREHOUSE_H

#include "Warehouse.h"
#include "MooseTypes.h"

#include <vector>
#include <map>
#include <string>
#include <list>
#include <set>

class ExactGeometry;

/**
 * Warehouse for storing ExactGeometries
 */
class ExactGeometryWarehouse : public Warehouse<ExactGeometry>
{
public:
  ExactGeometryWarehouse();
  virtual ~ExactGeometryWarehouse();

  // Setup /////
  void initialSetup();
  void timestepSetup();
  void residualSetup();
  void jacobianSetup();

  const std::vector<MooseSharedPointer<ExactGeometry> > & allExactGeometries() const { return _all_exact_geometries; }

  std::vector<MooseSharedPointer<ExactGeometry> > & activeBlockExactGeometries(SubdomainID block) { return _active_block_exact_geometries[block]; }

  /**
   * @param ExactGeometry being added
   */
  void addExactGeometry(MooseSharedPointer<ExactGeometry> & exact_geometry);

protected:
  /// all ExactGeometries
  std::vector<MooseSharedPointer<ExactGeometry> > _all_exact_geometries;

  /// nodal kernels active on a block
  std::map<SubdomainID, std::vector<MooseSharedPointer<ExactGeometry> > > _active_block_exact_geometries;
};

#endif // EXACTGEOMETRYWAREHOUSE_H
