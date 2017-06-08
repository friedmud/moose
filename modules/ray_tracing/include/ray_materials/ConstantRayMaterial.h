#ifndef CONSTANTRAYMATERIAL_H
#define CONSTANTRAYMATERIAL_H

// MOOSE Includes
#include "RayMaterial.h"

class ConstantRayMaterial;

template <>
InputParameters validParams<ConstantRayMaterial>();

/**
 * Holds all of the XS data for one material
 */
class ConstantRayMaterial : public RayMaterial
{
public:
  ConstantRayMaterial(const InputParameters & params);

  /// Sigma total for each energy group (starting with highest energy)
  const std::vector<double> & _sigma_t;
};

#endif
