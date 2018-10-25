//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef MESHSIDESETGENERATOR_H
#define MESHSIDESETGENERATOR_H

#include "MeshGenerator.h"

// Forward declarations
class MeshSideSetGenerator;

template <>
InputParameters validParams<MeshSideSetGenerator>();

/**
 * Add lower dimensional elements along the faces contained in a side set
 */
class MeshSideSetGenerator : public MeshGenerator
{
public:
  MeshSideSetGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate();

protected:
  std::unique_ptr<MeshBase> & _input;

  /// Block ID to assign to the region
  const subdomain_id_type _block_id;
};

#endif // MESHSIDESETGENERATOR_H
