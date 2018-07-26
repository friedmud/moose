#ifndef TRACERAY_H
#define TRACERAY_H

// MOOSE Includes
#include "Moose.h"
#include "MooseTypes.h"
#include "StaticallyAllocatedSet.h"
#include "PerfGraphInterface.h"

// libMesh Includes
#include "libmesh/id_types.h"

// System Includes
#include <map>

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

// Maximum number of point neighbors possible
// With normal topology the max is 8 for Hex8s
// Using 16 for safety
#define MAX_POINT_NEIGHBORS 16

extern const std::vector<std::vector<unsigned int>> quad4_side_to_children;

class TraceRay : public PerfGraphInterface
{
public:
  TraceRay(PerfGraph & perf_graph,
           const MeshBase & mesh,
           const BoundingBox & b_box,
           const std::map<const Elem *, std::vector<Point>> & elem_normals,
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

  /**
   * The number of intersections that were done
   */
  unsigned long intersections() { return _num_intersections; }

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
  const std::map<const Elem *, std::vector<Point>> & _elem_normals;
  unsigned int _halo_size;
  double _ray_max_distance;
  double _ray_length;

  /// Whether or not it's ok that a Ray fails tracing
  bool _tolerate_failure;

  THREAD_ID _tid;

  unsigned long _num_intersections;
  /*
    PerfID _trace_timer;
    PerfID _intersection_timer;
    PerfID _try_harder_timer;
    PerfID _found_intersection;
  */

  /// Which boundary IDs have already been applied, The '18' is the maximum.
  /// I expect that the most a ray will have is 6: 3 from one corner, 3 in another
  /// So I just multiplied by 3 to be extra safe!
  MooseUtils::StaticallyAllocatedSet<BoundaryID, 18> _applied_ids;

  /// For calling point_neighbors
  MooseUtils::StaticallyAllocatedSet<const Elem *, MAX_POINT_NEIGHBORS> _point_neighbors;

private:
  void find_point_neighbors(
      const Elem * current_elem,
      const Point & p,
      MooseUtils::StaticallyAllocatedSet<const Elem *, MAX_POINT_NEIGHBORS> & neighbor_set,
      const std::shared_ptr<Ray> &);

  void possiblyOnBoundary(const std::shared_ptr<Ray> & ray,
                          const Point & incoming_point,
                          const Elem * current_elem,
                          unsigned int incoming_side,
                          Point & intersection_point,
                          int & intersected_side);

  void tryToMoveThroughPointNeighbors(const Elem * current_elem,
                                      const ElemType elem_type,
                                      const Point & incoming_point,
                                      const std::shared_ptr<Ray> & ray,
                                      const Point & ray_direction,
                                      Point & intersection_point,
                                      Point & boundary_intersection_point,
                                      int & intersected_side,
                                      const Elem *& best_neighbor,
                                      unsigned int & best_side);

  void checkForCornerHitAndApplyBCs(const Elem * current_elem,
                                    const Point & intersection_point,
                                    std::shared_ptr<Ray> & ray);

  unsigned int _mesh_dim;

  Point _work_point;
  Point _work_point2;

  // Boundary IDs
  std::vector<BoundaryID> _ids;

  // Reusable
  std::vector<const Elem *> _active_neighbor_children;
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
