#ifndef TRACERAY_H
#define TRACERAY_H

// MOOSE Includes
#include "Moose.h"
#include "MooseTypes.h"

// libMesh Includes
#include "libmesh/id_types.h"

// forward declares
class Ray;
class MooseMesh;
class FEProblem;
class RayKernel;
class RayProblemBase;
class RaySystem;

class TraceRay
{
public:
  TraceRay(const MooseMesh & mesh,
           const BoundingBox & b_box,
           unsigned int halo_size,
           Real ray_max_distance,
           Real ray_length,
           THREAD_ID tid);

  void trace(std::shared_ptr<Ray> & ray);

protected:
  virtual void subdomainSetup(subdomain_id_type /* current_subdomain */,
                              std::shared_ptr<Ray> & /* ray */)
  {
  }

  virtual void onSegment(const Elem * /* current_elem */,
                         const Point & /* incoming_point */,
                         const Point & /* intersection_point */,
                         bool /* ends_in_elem */)
  {
  }

  virtual void onBoundary(boundary_id_type /* bid */,
                          const Elem * /* neighbor */,
                          unsigned int /* side */,
                          std::shared_ptr<Ray> & /* ray */)
  {
  }

  virtual void finishedBoundary() {}

  const MooseMesh & _mesh;
  BoundingBox _b_box;
  unsigned int _halo_size;
  double _ray_max_distance;
  double _ray_length;
  THREAD_ID _tid;

  const Elem * current_elem = NULL;
  int intersected_side = -1;
  Point intersection_point;
  Point boundary_intersection_point;
};

class RayProblemTraceRay : public TraceRay
{
public:
  RayProblemTraceRay(RayProblemBase & ray_problem,
                     const MooseMesh & mesh,
                     unsigned int halo_size,
                     Real ray_max_distance,
                     Real ray_length,
                     THREAD_ID tid);

  virtual void subdomainSetup(subdomain_id_type current_subdomain,
                              std::shared_ptr<Ray> & ray) override;

  virtual void onSegment(const Elem * current_elem,
                         const Point & incoming_point,
                         const Point & intersection_point,
                         bool ends_in_elem) override;

  virtual void onBoundary(boundary_id_type bid,
                          const Elem * neighbor,
                          unsigned int side,
                          std::shared_ptr<Ray> & ray) override;

  virtual void finishedBoundary() override;

protected:
  RayProblemBase & _ray_problem;
  RaySystem & _ray_system;

  const std::vector<MooseSharedPointer<RayKernel>> * _ray_kernels;

  std::set<BoundaryID> _applied_ids;
};

#endif
