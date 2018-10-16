//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CONCENTRICCIRCLEMESHGENERATOR_H
#define CONCENTRICCIRCLEMESHGENERATOR_H

#include "MeshGenerator.h"
#include "MooseEnum.h"

// Forward declarations
class ConcentricCircleMeshGenerator;

template <>
InputParameters validParams<ConcentricCircleMeshGenerator>();

/**
 * Generates an annular mesh given all the parameters
 */
class ConcentricCircleMeshGenerator : public MeshGenerator
{
public:
  ConcentricCircleMeshGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate();

protected:

  /// Number of sectors in one quadrant
  unsigned int _num_sectors;
  
  /// Radii of concentric circles
  std::vector<Real> _radii;
  
  /// Number of rings in each circle or in the moderator
  std::vector<unsigned int> _rings;
  
  /// Size of inner square in relation to radius of the innermost concentric circle
  Real _inner_mesh_fraction;
  
  /// Adding the moderator is optional
  bool _has_outer_square;
  Real _pitch;
  
  /// Volume preserving function is optional
  bool _preserve_volumes;
  
  /// Control of which portion of mesh will be developed
  MooseEnum _portion;
  
};

#endif // ANNULARMESHGENERATOR_H