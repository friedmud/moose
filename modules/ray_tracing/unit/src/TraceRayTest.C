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

// libMesh Includes

// Helper for comparing two vectors of doubles
#define ASSERT_VECTOR_DOUBLE_EQ(vec1, vec2, message)                                               \
  for (unsigned int vadei = beginIndex(vec1); vadei < vec1.size(); vadei++)                        \
    ASSERT_NEAR(vec1[vadei], vec2[vadei], TOLERANCE) << message;

class DebugTraceRay : public TraceRay
{
public:
  DebugTraceRay(const MeshBase & mesh,
                const BoundingBox & b_box,
                unsigned int halo_size,
                Real ray_max_distance,
                Real ray_length,
                THREAD_ID tid,
                bool verbose = false)
    : TraceRay(mesh, b_box, halo_size, ray_max_distance, ray_length, tid), _verbose(verbose)
  {
  }

  std::vector<dof_id_type> _elements_traversed;
  std::vector<boundary_id_type> _boundaries_hit;
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

    if (_verbose)
    {
      std::cout << "Current Elem: " << current_elem->id() << std::endl;
      std::cout << "Segment length: " << segment_length << std::endl;
    }
  }

  virtual void onBoundary(const Elem * /* current_elem */,
                          const Point & /* intersection_point */,
                          unsigned int /* intersected_side */,
                          boundary_id_type bid,
                          const Elem * /* neighbor */,
                          std::shared_ptr<Ray> & ray)
  {
    if (_verbose)
      std::cout << "Hit boundary " << bid << std::endl;

    _boundaries_hit.push_back(bid);

    ray->setShouldContinue(false);
  }

  virtual void finishedBoundary() {}

protected:
  bool _verbose;
};

TEST_F(TraceRayTest, trace_single_2d_quad)
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

TEST_F(TraceRayTest, trace_adapted_2d_quads)
{
  InputParameters mesh_params = _factory->getValidParams("GeneratedMesh");
  mesh_params.set<MooseEnum>("dim") = "2";
  mesh_params.set<std::string>("_object_name") = "mesh";

  mesh_params.set<Real>("xmin") = 0.0;
  mesh_params.set<Real>("ymin") = 0.0;

  mesh_params.set<Real>("xmax") = 2.0;
  mesh_params.set<Real>("ymax") = 1.0;

  mesh_params.set<unsigned int>("nx") = 2;
  mesh_params.set<unsigned int>("ny") = 1;

  std::unique_ptr<MooseMesh> mesh = libmesh_make_unique<GeneratedMesh>(mesh_params);
  mesh->buildMesh();

  // Setup refinement pattern
  mesh->elemPtr(0)->set_refinement_flag(Elem::DO_NOTHING);
  mesh->elemPtr(1)->set_refinement_flag(Elem::REFINE);

  MeshRefinement refinement(*mesh);
  refinement.refine_and_coarsen_elements();

  // Trace a Ray from left side to right.  Should go through two of the new children.
  auto ray = std::make_shared<Ray>(Point(0, 0.2, 0), Point(2, 0.2, 0), 1, mesh->elemPtr(0), 3);

  ray->setID(0);

  DebugTraceRay tr(*mesh, BoundingBox(Point(0, 0, 0), Point(2, 1, 0)), 1, 1e5, 1e5, 0);

  tr.trace(ray);

  ASSERT_TRUE(tr._elements_traversed == std::vector<dof_id_type>({0, 2, 3}))
      << "Correct elements not traversed!";

  ASSERT_VECTOR_DOUBLE_EQ(
      tr._segment_lengths, std::vector<Real>({1.0, 0.5, 0.5}), "Segment lengths don't match!")
}

TEST_F(TraceRayTest, trace_more_adapted_2d_quads)
{
  InputParameters mesh_params = _factory->getValidParams("GeneratedMesh");
  mesh_params.set<MooseEnum>("dim") = "2";
  mesh_params.set<std::string>("_object_name") = "mesh";

  mesh_params.set<Real>("xmin") = 0.0;
  mesh_params.set<Real>("ymin") = 0.0;

  mesh_params.set<Real>("xmax") = 2.0;
  mesh_params.set<Real>("ymax") = 1.0;

  mesh_params.set<unsigned int>("nx") = 2;
  mesh_params.set<unsigned int>("ny") = 1;

  std::unique_ptr<MooseMesh> mesh = libmesh_make_unique<GeneratedMesh>(mesh_params);
  mesh->buildMesh();

  // Setup refinement pattern
  MeshRefinement refinement(*mesh);
  refinement.clean_refinement_flags();

  mesh->elemPtr(1)->set_refinement_flag(Elem::REFINE);
  refinement.refine_and_coarsen_elements();

  refinement.clean_refinement_flags();

  // Refinement 2
  mesh->elemPtr(2)->set_refinement_flag(Elem::REFINE);
  refinement.refine_and_coarsen_elements();

  // Trace a Ray from left side to right.  Should go through two of the new children.
  auto ray = std::make_shared<Ray>(Point(0, 0.2, 0), Point(2, 0.2, 0), 1, mesh->elemPtr(6), 3);

  ray->setID(0);

  DebugTraceRay tr(*mesh, BoundingBox(Point(0, 0, 0), Point(2, 1, 0)), 1, 1e5, 1e5, 0);
  /*
    for (auto elem_it = mesh->getMesh().elements_begin(); elem_it != mesh->getMesh().elements_end();
         elem_it++)
      std::cout << *(*elem_it) << std::endl;
  */

  tr.trace(ray);

  ASSERT_TRUE(tr._elements_traversed == std::vector<dof_id_type>({6, 7, 10, 11, 3}))
      << "Correct elements not traversed!";

  ASSERT_VECTOR_DOUBLE_EQ(tr._segment_lengths,
                          std::vector<Real>({0.5, 0.5, 0.25, 0.25, 0.5}),
                          "Segment lengths don't match!")
}

TEST_F(TraceRayTest, trace_adapted_2d_quads_strike_level_mismatch_corner)
{
  InputParameters mesh_params = _factory->getValidParams("GeneratedMesh");
  mesh_params.set<MooseEnum>("dim") = "2";
  mesh_params.set<std::string>("_object_name") = "mesh";

  mesh_params.set<Real>("xmin") = 0.0;
  mesh_params.set<Real>("ymin") = 0.0;

  mesh_params.set<Real>("xmax") = 2.0;
  mesh_params.set<Real>("ymax") = 1.0;

  mesh_params.set<unsigned int>("nx") = 2;
  mesh_params.set<unsigned int>("ny") = 1;

  std::unique_ptr<MooseMesh> mesh = libmesh_make_unique<GeneratedMesh>(mesh_params);
  mesh->buildMesh();

  // Setup refinement pattern
  mesh->elemPtr(0)->set_refinement_flag(Elem::DO_NOTHING);
  mesh->elemPtr(1)->set_refinement_flag(Elem::REFINE);

  MeshRefinement refinement(*mesh);
  refinement.refine_and_coarsen_elements();

  // Trace a Ray from left side to right.  Should go through two of the new children.
  auto ray = std::make_shared<Ray>(Point(0, 0, 0), Point(2, 1, 0), 1, mesh->elemPtr(0), 3);

  ray->setID(0);

  DebugTraceRay tr(*mesh, BoundingBox(Point(0, 0, 0), Point(2, 1, 0)), 1, 1e5, 1e5, 0, false);

  tr.trace(ray);

  ASSERT_TRUE(tr._elements_traversed == std::vector<dof_id_type>({0, 4, 5}))
      << "Correct elements not traversed!";

  ASSERT_TRUE(tr._boundaries_hit == std::vector<boundary_id_type>({1, 2}))
      << "Correct boundaries not hit!" << std::endl;

  ASSERT_VECTOR_DOUBLE_EQ(tr._segment_lengths,
                          std::vector<Real>({1.1180339887, 1.1180339887 / 2., 1.1180339887 / 2.}),
                          "Segment lengths don't match!")
}
