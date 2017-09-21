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
#include "TraceRayTest.h"

// RayTracing Includes
#include "Ray.h"
#include "TraceRay.h"

class DebugTraceRay : public TraceRay
{
public:
  DebugTraceRay(const MeshBase & mesh,
                const BoundingBox & b_box,
                unsigned int halo_size,
                Real ray_max_distance,
                Real ray_length,
                THREAD_ID tid)
    : TraceRay(mesh, b_box, halo_size, ray_max_distance, ray_length, tid)
  {
  }

  std::vector<dof_id_type> _elements_traversed;
  std::vector<Real> _segment_lengths;

protected:
  virtual void subdomainSetup(subdomain_id_type /* current_subdomain */,
                              std::shared_ptr<Ray> & /* ray */)
  {
  }

  virtual void onSegment(const Elem * current_elem,
                         const Point & incoming_point,
                         const Point & intersection_point,
                         bool /* ends_in_elem */)
  {
    Real segment_length = (intersection_point - incoming_point).size();

    _elements_traversed.push_back(current_elem->id());
    _segment_lengths.push_back(segment_length);

    std::cout << "Current Elem: " << current_elem->id() << std::endl;
    std::cout << "Segment length: " << segment_length << std::endl;
  }

  virtual void onBoundary(const Elem * /* current_elem */,
                          const Point & /* intersection_point */,
                          unsigned int /* intersected_side */,
                          boundary_id_type /* bid */,
                          const Elem * /* neighbor */,
                          std::shared_ptr<Ray> & ray)
  {
    std::cout << "Hit boundary" << std::endl;

    ray->setShouldContinue(false);
  }

  virtual void finishedBoundary() {}
};

// Helper for comparing two vectors of doubles
#define ASSERT_VECTOR_DOUBLE_EQ(vec1, vec2, message)                                               \
  for (unsigned int vadei = beginIndex(vec1); vadei < vec1.size(); vadei++)                        \
    ASSERT_DOUBLE_EQ(vec1[vadei], vec2[vadei]) << message;

TEST_F(TraceRayTest, trace_single)
{
  InputParameters mesh_params = _factory->getValidParams("GeneratedMesh");
  mesh_params.set<MooseEnum>("dim") = "2";
  mesh_params.set<std::string>("_object_name") = "mesh";

  mesh_params.set<Real>("xmin") = 0.0;
  mesh_params.set<Real>("ymin") = 0.0;

  mesh_params.set<Real>("xmax") = 1.0;
  mesh_params.set<Real>("ymax") = 1.0;

  mesh_params.set<unsigned int>("nx") = 1;
  mesh_params.set<unsigned int>("ny") = 1;

  std::unique_ptr<MooseMesh> mesh = libmesh_make_unique<GeneratedMesh>(mesh_params);
  mesh->buildMesh();

  auto ray = std::make_shared<Ray>(Point(0, 0.5, 0), Point(1, 0.5, 0), 1, mesh->elemPtr(0), 3);

  ray->setID(0);

  DebugTraceRay tr(*mesh, BoundingBox(Point(0, 0, 0), Point(1, 1, 1)), 1, 1e5, 1e5, 0);

  tr.trace(ray);

  ASSERT_TRUE(tr._elements_traversed == std::vector<dof_id_type>({0}))
      << "Element 0 not traversed!";

  ASSERT_VECTOR_DOUBLE_EQ(
      tr._segment_lengths, std::vector<Real>({1.0}), "Segment lengths don't match!")
}
