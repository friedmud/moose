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
#include "RayTest.h"

#include "Ray.h"

TEST_F(RayTest, pack_unpack)
{
  auto original =
      std::make_shared<Ray>(Point(0.1, 0.2, 0.3), Point(0.4, 0.5, 0.6), 4, _mesh->elemPtr(0), 1);

  original->_data[0] = 1.2;
  original->_data[1] = 1.3;
  original->_data[2] = 1.4;
  original->_data[3] = 1.5;

  original->processorCrossings() = 14;
  original->intersections() = 17;
  original->distance() = 2.7;
  original->integratedDistance() = 2.1;

  original->setAzimuthalSpacing(0.03);
  original->setAzimuthalWeight(0.4);
  original->setPolarSpacing(0.5);
  original->setPolarSins({0.6, 0.7});
  original->setPolarWeights({0.8, 0.9});

  std::vector<Real> packed_data;

  std::back_insert_iterator<std::vector<Real>> back_inserter(packed_data);

  libMesh::Parallel::Packing<std::shared_ptr<Ray>>::pack(original, back_inserter, NULL);

  auto unpacked = libMesh::Parallel::Packing<std::shared_ptr<Ray>>::unpack(packed_data.begin(),
                                                                           &_mesh->getMesh());

  ASSERT_TRUE((*original) == (*unpacked)) << "Original Ray differs from the packed Ray";
}
