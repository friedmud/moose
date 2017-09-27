#ifndef TRACERAY_H
#define TRACERAY_H

// MOOSE Includes
#include "Moose.h"
#include "MooseTypes.h"

// libMesh Includes
#include "libmesh/id_types.h"

// forward declares
class Ray;
class FEProblem;
class RayKernel;
class RayProblemBase;
class RaySystem;

namespace libMesh
{
class Mesh;
}

extern const std::vector<std::vector<unsigned int>> quad4_side_to_children;

class TraceRay
{
public:
  TraceRay(const MeshBase & mesh,
           const BoundingBox & b_box,
           unsigned int halo_size,
           Real ray_max_distance,
           Real ray_length,
           bool tolerate_failure,
           THREAD_ID tid);

  void trace(std::shared_ptr<Ray> & ray);

  /**
   * Find the child of a refined element that contains the given point
   *
   * @return The element containing the point
   *
   * @param current_child The children of this element will be searched
   * @param intersection_point The point to be searched for
   * @param children_ids The embedding matrix ID of the childred that should be searched.  Useful
   * for narrowing down the correct children to search when you can rule some out (for instance: if
   * you only need to search the children on one specific side of the parent element.  Be aware:
   * these are used all the way down - so make sure what you do makes sense!
   */
  static Elem * recursivelyFindChildContainingPoint(const Elem * current_child,
                                                    const Point & intersection_point,
                                                    const std::vector<unsigned int> & children_ids);

  /**
   * Find the child that contains the point on a side of the current element
   *
   * This is an optimization over just using recursivelyFindCHildContinaingPoint that only looks in
   * children that are on the side specified.
   *
   * @param current_elem The elem to search
   * @param p The physical point to look for
   * @param side The side of current_elem that should be searched
   */
  static Elem * childOnSide(const Elem * current_elem, const Point & p, unsigned int side);

  /**
   * Get the neighbor on the intersected_side of current_elem.
   *
   * In the case of mesh adaptivity this will return the most refined element that is on the side of
   * current_elem and contains the intersection_point
   *
   * @param current_elem The elem you want to find the neighbor for
   * @param intersected_side The side of current_elem you want a neighbor for
   * @param intersection_point The point on the side you want the neighbor for (only used when mesh
   * adaptivity is active to find the child that contains the point)
   */
  static Elem *
  getNeighbor(const Elem * current_elem, unsigned int intersected_side, Point & intersection_point);

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

  virtual void onBoundary(const Elem * /* current_elem */,
                          const Point & /* intersection_point */,
                          unsigned int /* intersected_side */,
                          boundary_id_type /* bid */,
                          const Elem * /* neighbor */,
                          std::shared_ptr<Ray> & /* ray */)
  {
  }

  virtual void finishedBoundary() {}

  const MeshBase & _mesh;
  BoundingBox _b_box;
  unsigned int _halo_size;
  double _ray_max_distance;
  double _ray_length;

  /// Whether or not it's ok that a Ray fails tracing
  bool _tolerate_failure;

  THREAD_ID _tid;

  /// Which boundary IDs have already been applied
  std::set<BoundaryID> _applied_ids;
};

class RayProblemTraceRay : public TraceRay
{
public:
  RayProblemTraceRay(RayProblemBase & ray_problem,
                     const MeshBase & mesh,
                     unsigned int halo_size,
                     Real ray_max_distance,
                     Real ray_length,
                     bool tolerate_failure,
                     THREAD_ID tid);

  virtual void subdomainSetup(subdomain_id_type current_subdomain,
                              std::shared_ptr<Ray> & ray) override;

  virtual void onSegment(const Elem * current_elem,
                         const Point & incoming_point,
                         const Point & intersection_point,
                         bool ends_in_elem) override;

  virtual void onBoundary(const Elem * current_elem,
                          const Point & intersection_point,
                          unsigned int intersected_side,
                          boundary_id_type bid,
                          const Elem * /* neighbor */,
                          std::shared_ptr<Ray> & ray) override;

  virtual void finishedBoundary() override;

protected:
  RayProblemBase & _ray_problem;
  RaySystem & _ray_system;

  const std::vector<MooseSharedPointer<RayKernel>> * _ray_kernels;
};

#endif
