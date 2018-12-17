#include "TraceRay.h"

// Local Includes
#include "RayKernel.h"
#include "RayBC.h"
#include "RayProblem.h"

// Moose Includes
#include "RayTracing.h"
#include "Ray.h"
#include "MooseError.h"
#include "FEProblem.h"
#include "RayKernel.h"

// libMesh Includes
#include "libmesh/plane.h"
#include "libmesh/point.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"
#include "libmesh/remote_elem.h"

// System Includes
#include <list>
#include <unistd.h>

unsigned long int ray_count = 0;

// 0 ray_count: 156007 ray_id: 7244 Endless loop!

unsigned long int debug_ray = 216;

unsigned long int debug_ray_id = 42846600;

unsigned long int debug_ray_pid = 1;

//#define DEBUG_IF ray_count == debug_ray

#define DEBUG_IF ray->id() == debug_ray_id // ray_count > 197000 && debug_ray_id == ray->id()

// && ray_count == debug_ray

// #define USE_DEBUG_RAY

const std::vector<std::vector<unsigned long int>> quad4_side_to_children = {
    {0, 1}, {1, 3}, {3, 2}, {2, 0}};

// The sides that adjoin each edge in a hex8
const std::vector<std::vector<unsigned int>> hex8_edge_to_sides =
{
//0      1       2       3       4       5       6       7       8       9       10      11
  {0,1}, {0, 2}, {0, 3}, {0, 4}, {1, 4}, {1, 2}, {2, 3}, {3, 4}, {1, 5}, {2, 5}, {3, 5}, {4, 5}
};







TraceRay::TraceRay(PerfGraph & perf_graph,
                   const MeshBase & mesh,
                   const BoundingBox & b_box,
                   const std::map<const Elem *, std::vector<Point>> & elem_normals,
                   unsigned long int halo_size,
                   Real ray_max_distance,
                   Real ray_length,
                   bool tolerate_failure,
                   THREAD_ID tid)
  : PerfGraphInterface(perf_graph, "TraceRay"),
    _mesh(mesh),
    _b_box(b_box),
    _elem_normals(elem_normals),
    _halo_size(halo_size),
    _ray_max_distance(ray_max_distance),
    _ray_length(ray_length),
    _tolerate_failure(tolerate_failure),
    _tid(tid),
    _num_intersections(0)
/*
_trace_timer(registerTimedSection("trace", 1)),
_intersection_timer(registerTimedSection("intersection", 1)),
_try_harder_timer(registerTimedSection("try_harder", 1)),
_found_intersection(registerTimedSection("found_intersection", 1))
*/
{
  _ids.reserve(20);
}

RayProblemTraceRay::RayProblemTraceRay(RayProblemBase & ray_problem,
                                       const MeshBase & mesh,
                                       unsigned long int halo_size,
                                       Real ray_max_distance,
                                       Real ray_length,
                                       bool tolerate_failure,
                                       THREAD_ID tid)
  : TraceRay(ray_problem.perfGraph(),
             mesh,
             ray_problem.boundingBox(),
             ray_problem.elemNormals(),
             halo_size,
             ray_max_distance,
             ray_length,
             tolerate_failure,
             tid),
    _ray_problem(ray_problem),
    _ray_system(ray_problem.raySystem())
{
}

void
RayProblemTraceRay::subdomainSetup(subdomain_id_type current_subdomain, std::shared_ptr<Ray> & ray)
{
  _ray_system.subdomainSetup(current_subdomain, _tid);
  _ray_kernels = &_ray_system.getRayKernels(current_subdomain, _tid);

  for (auto & ray_kernel : *_ray_kernels)
    ray_kernel->setRay(ray);
}

void
RayProblemTraceRay::onSegment(const Elem * current_elem,
                              const Point & incoming_point,
                              const Point & intersection_point,
                              bool ends_in_elem)
{
  // The "true" causes only sigmaT to be reinited
  _ray_system.reinitElem(current_elem, _tid, true);

  for (auto & ray_kernel : *_ray_kernels)
    ray_kernel->onSegment(current_elem, incoming_point, intersection_point, ends_in_elem);
}

void
RayProblemTraceRay::onBoundary(const Elem * current_elem,
                               const Point & intersection_point,
                               unsigned long int intersected_side,
                               boundary_id_type bid,
                               const Elem * /* neighbor */,
                               std::shared_ptr<Ray> & ray)
{
  if (_ray_system.hasRayBCs(bid, _tid))
  {
    auto & ray_bcs = _ray_system.getRayBCs(bid, _tid);

    for (auto & bc : ray_bcs)
    {
#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        std::cerr << ray_count << " Applying BC!" << std::endl;
#endif

      bc->apply(current_elem, intersected_side, intersection_point, ray);
    }
  }
}

void
RayProblemTraceRay::finishedBoundary()
{}

/**
 * Pass in:
 * Ray: o -> d
 * Side: v0 -> v1
 *
 * Internally:
 * Ray: p -> t*r
 * Side: q -> u*s
 *
 * From: https://stackoverflow.com/a/565282
 *
 * Note this is the 3D version of the 2D only code below
 *
 * Also taken from "Intersection of two lines in three-space" by Ronald Goldman in Graphics Gems on page 304
 */
bool
lineLineIntersect3DVanilla(const std::shared_ptr<Ray> &

#ifdef USE_DEBUG_RAY
                           ray
#endif
                           ,
                           const Point & o,
                           const Point & d,
                           const Point & v0,
                           const Point & v1,
                           Point & intersection_point,
                           Real & intersection_distance)
{
  const Real EPSILON = 1e-8;

  const Point & p = o;
  const Point & q = v0;

  auto r = (d - o);
  auto r_size = r.size();
  r /= r_size;

  auto s = (v1 - v0);
  auto s_size = s.size();
  s /= s_size;;

  auto rxs = r.cross(s);

  auto rxs_size_sq = rxs.size_sq();

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << " rxs_size_sq: " << rxs_size_sq << std::endl;
#endif

  if (rxs_size_sq < 1e-10) // Lines are parallel or colinear
    return false;

  auto qmp = q - p;

  // Need to make two TypeTensors here to take their determinant
  // Load into columns: {qmp, s, rxs}
  TypeTensor<Real> detfinder;
  detfinder.slice(0) = qmp;
  detfinder.slice(1) = s;
  detfinder.slice(2) = rxs;

  auto t = detfinder.det() / rxs_size_sq;

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
      std::cerr << " t: " << t << std::endl;
#endif

  if (t < -EPSILON || t > r_size + EPSILON)
    return false;

  // Change the second column to be r
  detfinder.slice(1) = r;

  auto u = detfinder.det() / rxs_size_sq;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << " u: " << u << std::endl;
#endif

  if (u < -EPSILON || u > s_size + EPSILON)
    return false;

  Point ppoint = p + (t * r);
  Point qpoint = q + (u * s);

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
    {
      std::cerr << " ppoint: " << ppoint << std::endl;
      std::cerr << " qpoint: " << qpoint << std::endl;
      std::cerr << " difference: " << ppoint - qpoint << std::endl;
    }
#endif

  if (ppoint.absolute_fuzzy_equals(qpoint, 1e-8))
  {
    intersection_point = p + (t * r);
    intersection_distance = t;
    return true;
  }

  // Not parallel, but don't intersect
  return false;
}

/**
 * Pass in:
 * Ray: o -> d
 * Side: v0 -> v1
 *
 * Internally:
 * Ray: p -> t*r
 * Side: q -> u*s
 *
 * From: https://stackoverflow.com/a/565282
 */
bool
lineLineIntersect2DVanilla(const Point & o,
                           const Point & d,
                           const Point & v0,
                           const Point & v1,
                           Point & intersection_point,
                           Real & intersection_distance)

{
  const Point & p = o;
  const Point & q = v0;

  auto r = (d - o);
  auto s = (v1 - v0);

  auto rxs = r(0) * s(1) - r(1) * s(0);

  if (std::abs(rxs) < 1e-10) // Lines are parallel or colinear
    return false;

  auto qmp = q - p;

  auto t = (qmp(0) * s(1) - qmp(1) * s(0)) / rxs;
  auto u = (qmp(0) * r(1) - qmp(1) * r(0)) / rxs;

  if ((0 < t + 4e-9 && t - 4e-9 <= 1.0) && (0 < u + 4e-9 && u - 4e-9 <= 1.0)) // Do they intersect
  {
    intersection_point = p + (t * r);
    intersection_distance = t;
    return true;
  }

  // Not parallel, but don't intersect
  return false;
}

/**
 * Derived from: http://stackoverflow.com/a/565282/2042320
 * q -> q + u*s: Ray
 * p -> p + t*r: Side
 */
template <typename T>
int
sideIntersectedByLine2D(const Elem * current_elem,
                        unsigned long int incoming_side,
                        const Point & incoming_point,
                        const std::shared_ptr<Ray> & ray,
                        const Point & ray_direction,
                        const std::map<const Elem *, std::vector<Point>> & /*elem_normals*/,
                        Point & intersection_point,
                        Point & /*boundaryintersection_point*/)
{
  int intersected_side = -1;
  unsigned long int n_sides = current_elem->n_sides();

  Point current_intersection_point;
  double current_intersection_distance = 0;

  {
    //    int boundary_side = -1;

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
    {
      std::cerr << "Incoming point: " << incoming_point << std::endl;
      std::cerr << "Incoming side: " << incoming_side << std::endl;
      std::cerr << "Ray end: " << ray->end() << std::endl;
    }

#endif

    // Bump the starting point down the path a bit
    Point starting_point = incoming_point + (1e-9 * ray_direction);

    double intersection_distance = 0;

    // Unfortunately - backface culling ended up being slower than just testing the intersection
    // *shrug*
    //    auto & normals = elem_normals.at(current_elem);

    for (unsigned long int i = 0; i < n_sides; i++)
    {
      if (i == incoming_side) // Don't search backwards
        continue;

        // Backface culling
        //      if (normals[i] * ray_direction < -TOLERANCE)
        //        continue;

#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
      {
        std::cerr << "Side: " << i << std::endl;
        std::cerr << " Side nodes: " << std::endl;
        std::cerr << "  " << current_elem->point(T::side_nodes_map[i][0]) << std::endl;
        std::cerr << "  " << current_elem->point(T::side_nodes_map[i][1]) << std::endl;
      }
#endif

      auto intersected = lineLineIntersect2DVanilla(starting_point,
                                                    ray->end(),
                                                    current_elem->point(T::side_nodes_map[i][0]),
                                                    current_elem->point(T::side_nodes_map[i][1]),
                                                    current_intersection_point,
                                                    current_intersection_distance);

      // Do they intersect and is it further down the path than any other intersection
      if (intersected && current_intersection_distance > intersection_distance)
      {
        intersected_side = i;
        intersection_point = current_intersection_point;
        intersection_distance = current_intersection_distance;
      }
    }
  }

  return intersected_side;
}

// Code adapted from:
// http://webserver2.tecgraf.puc-rio.br/~mgattass/cg/trbRR/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf
bool
rayIntersectsTriangleCulled(const Point & O,
                            const Point & D,
                            const Point & V00,
                            const Point & V10,
                            const Point & V11,
                            Real & u,
                            Real & v,
                            Real & t,
                            const std::shared_ptr<Ray> &
#ifdef USE_DEBUG_RAY
                                ray
#endif
)
{
  const Real EPSILON = 0.000001;
  const Point & vert0 = V00;
  const Point & vert1 = V10;
  const Point & vert2 = V11;

  const Point & orig = O;
  const Point & dir = D;

  auto edge1 = vert1 - vert0;
  auto edge2 = vert2 - vert0;

  auto pvec = dir.cross(edge2);

  auto det = edge1 * pvec;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "  det: " << det << std::endl;
#endif

  if (det < EPSILON)
    return false;

  auto tvec = orig - vert0;

  u = tvec * pvec;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "  u: " << u << std::endl;
#endif

  if (u < 0. || u > det)
    return false;

  auto qvec = tvec.cross(edge1);

  v = dir * qvec;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "  v: " << v << std::endl;
#endif

  if (v < 0. || u + v > det)
    return false;

  t = edge2 * qvec;

  auto inv_det = 1. / det;

  t *= inv_det;
  u *= inv_det;
  v *= inv_det;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "  t: " << v << std::endl;
#endif

  if (t <= 0)
    return false;

  return true;
}

// Code from: https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
bool
rayIntersectsTriangle(const Point & O,
                      const Point & D,
                      const Point & V00,
                      const Point & V10,
                      const Point & V11,
                      Real & u,
                      Real & v,
                      Real & t,
                      const std::shared_ptr<Ray> &
#ifdef USE_DEBUG_RAY
                          ray
#endif
                      /*
  Vector3D rayOrigin,
                           Vector3D rayVector,
                           Triangle* inTriangle,
                           Vector3D& outIntersectionPoint*/)
{
  const Real EPSILON = 0.000001;
  const Point & vertex0 = V00;
  const Point & vertex1 = V10;
  const Point & vertex2 = V11;

  const Point & rayOrigin = O;
  const Point & rayVector = D;

  Point edge1, edge2, h, s, q;

  Real a, f;

  edge1 = vertex1 - vertex0;
  edge2 = vertex2 - vertex0;
  h = rayVector.cross(edge2);
  a = edge1 * h;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "  a: " << a << std::endl;
#endif

  if (a > -EPSILON && a < EPSILON)
    return false;

  f = 1 / a;
  s = rayOrigin - vertex0;
  u = f * (s * h);

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "  u: " << u << std::endl;
#endif

  if (u < -EPSILON || u > 1. + EPSILON)
    return false;
  q = s.cross(edge1);

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "  q: " << q << std::endl;
#endif

  v = f * rayVector * q;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "  v: " << v << std::endl;
#endif

  if (v < -EPSILON || u + v > 1. + EPSILON)
    return false;

  // At this stage we can compute t to find out where the intersection point is on the line.
  t = f * edge2 * q;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "  t: " << t << std::endl;
#endif

  if (t > 0) // ray intersection
    return true;
  else // This means that there is a line intersection but not a ray intersection.
    return false;
}

bool
intersectQuadUsingTriangles(const Point & O,
                            const Point & D,
                            const Point & V00,
                            const Point & V10,
                            const Point & V11,
                            const Point & V01,
                            Real & u,
                            Real & v,
                            Real & t,
                            const std::shared_ptr<Ray> & ray)
{
  /*
  auto normal = (V10 - V00).cross(V11 - V10);

  // Backface culling
  if (D * normal > -TOLERANCE)
    return false;
  */

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << " Trying first triangle " << std::endl;
#endif
  if (rayIntersectsTriangleCulled(O, D, V00, V10, V11, u, v, t, ray))
    return true;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << " Trying second triangle " << std::endl;
#endif
  if (rayIntersectsTriangleCulled(O, D, V11, V01, V00, u, v, t, ray))
    return true;

  /*
#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << " Trying third triangle " << std::endl;
#endif
  if (rayIntersectsTriangle(O, D, V10, V11, V01, u, v, t, ray))
    return true;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << " Trying fourth triangle " << std::endl;
#endif
  if (rayIntersectsTriangle(O, D, V01, V00, V10, u, v, t, ray))
    return true;
  */

  return false;
}

// https://people.cs.kuleuven.be/~ares.lagae/publications/LD05ERQIT/LD05ERQIT_code.cpp
// An Efficient Ray-Quadrilateral Intersection Test
// Ares Lagae Philip Dutr ́e
// Department of Computer Science Katholieke Universiteit Leuven
// http://graphics.cs.kuleuven.be/publications/LD05ERQIT/LD05ERQIT_paper.pdf
bool
intersectQuad(const Point & O,
              const Point & D,
              const Point & V00,
              const Point & V10,
              const Point & V11,
              const Point & V01,
              Real & u,
              Real & v,
              Real & t,
              const std::shared_ptr<Ray> &
#ifdef USE_DEBUG_RAY
                  ray
#else
/* ray */
#endif
)
{
  // Reject rays using the barycentric coordinates of // the intersection point with respect to T.
  auto E01 = V10;
  E01 -= V00;

  auto E03 = V01;
  E03 -= V00;

  auto P = D.cross(E03);

  auto det = E01 * P;

  //  libMesh::err<<"Count: "<<count<<std::endl;

  if (std::abs(det) < TOLERANCE)
  {
#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
      libMesh::err << "1 Rejecting because " << std::abs(det) << " < " << TOLERANCE << std::endl;
#endif
    return false;
  }

  auto T = O;
  T -= V00;

  auto alpha = (T * P) / det;

  if (alpha < -1e-12)
  {
#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
      libMesh::err << "2 Rejecting because " << alpha << " < " << 0 << std::endl;
#endif
    return false;
  }

  /*
  if (alpha > 1)
  {

      libMesh::err<<"3 Rejecting because "<<alpha<<" > "<<1<<std::endl;

    return false;
  }
  */

  auto Q = T.cross(E01);

  // Compute the ray parameter of the intersection // point.
  t = (E03 * Q) / det;

  if (t < -1e-12)
  {
#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
      libMesh::err << "3 Rejecting because " << t << " < " << 0 << std::endl;
#endif

    return false;
  }

  auto beta = (D * Q) / det;

  if (beta < -1e-12)
  {
#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
      libMesh::err << "4 Rejecting because " << beta << " < " << 0 << std::endl;
#endif

    return false;
  }

  /*
    if (beta > 1)
    {

        libMesh::err<<"5 Rejecting because "<<beta<<" > "<<1+TOLERANCE<<std::endl;

      return false;
    }
  */

  // Reject rays using the barycentric coordinates of // the intersection point with respect to
  // T′.
  if ((alpha + beta) > 1)
  {
    auto E23 = V01;
    E23 -= V11;

    auto E21 = V10;
    E21 -= V11;

    auto P_prime = D.cross(E21);

    auto det_prime = E23 * P_prime;

    if (std::abs(det_prime) < TOLERANCE)
    {
#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        libMesh::err << "5 Rejecting because " << std::abs(det_prime) << " < " << TOLERANCE
                     << std::endl;
#endif

      return false;
    }

    auto T_prime = O;
    T_prime -= V11;

    auto alpha_prime = (T_prime * P_prime) / det_prime;

    if (alpha_prime < -1e-12)
    {
#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        libMesh::err << "6 Rejecting because " << alpha_prime << " < " << 0 << std::endl;
#endif

      return false;
    }

    auto Q_prime = T_prime.cross(E23);

    auto beta_prime = (D * Q_prime) / det_prime;

    if (beta_prime < -1e-12)
    {
#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        libMesh::err << "7 Rejecting because " << beta_prime << " < " << 0 << std::endl;
#endif

      return false;
    }
  }

  // Compute the barycentric coordinates of V11. E02 = V11 - V00
  auto E02 = V11;
  E02 -= V00;

  auto N = E01.cross(E03);

  Real alpha11;
  Real beta11;

  if ((std::abs(N(0)) >= std::abs(N(1))) && (std::abs(N(0)) >= std::abs(N(2))))
  {
    alpha11 = (E02(1) * E03(2) - E02(2) * E03(1)) / N(0);
    beta11 = (E01(1) * E02(2) - E01(2) * E02(1)) / N(0);
  }
  else if ((std::abs(N(1)) >= std::abs(N(0))) && (std::abs(N(1)) >= std::abs(N(2))))
  {
    alpha11 = (E02(2) * E03(0) - E02(0) * E03(2)) / N(1);
    beta11 = (E01(2) * E02(0) - E01(0) * E02(2)) / N(1);
  }
  else
  {
    alpha11 = (E02(0) * E03(1) - E02(1) * E03(0)) / N(2);
    beta11 = (E01(0) * E02(1) - E01(1) * E02(0)) / N(2);
  }

  // Compute the bilinear coordinates of the // intersection point.
  if (std::abs(alpha11 - 1) < TOLERANCE)
  {
    u = alpha;
    if (std::abs(beta11 - 1) < TOLERANCE)
      v = beta;
    else
      v = beta / (u * (beta11 - 1) + 1);
  }
  else if (std::abs(beta11 - 1) < TOLERANCE)
  {
    v = beta;
    u = alpha / (v * (alpha11 - 1) + 1);
  }
  else
  {
    auto A = -(beta11 - 1);
    auto B = alpha * (beta11 - 1) - beta * (alpha11 - 1) - 1;
    auto C = alpha;

    auto delta = (B * B) - (4 * A * C);
    auto Q = -0.5 * (B + ((B < 0.0 ? -1.0 : 1.0) * std::sqrt(delta)));
    u = Q / A;
    if ((u < 0) || (u > 1))
      u = C / Q;
    v = beta / (u * (beta11 - 1) + 1);
  }
  return true;
}

int
sideIntersectedByLineHex8(const Elem * current_elem,
                          unsigned long int incoming_side,
                          const Point & incoming_point,
                          const std::shared_ptr<Ray> & ray,
                          const Point & ray_direction,
                          const std::map<const Elem *, std::vector<Point>> & /*elem_normals*/,
                          Point & intersection_point,
                          Point & /*boundaryintersection_point*/)
{
  int intersected_side = -1;
  unsigned long int n_sides = 6;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "Incoming point: " << incoming_point << std::endl;
#endif

  {
    //    int boundary_side = -1;

    /*
    double q0 = incoming_point(0);
    double q1 = incoming_point(1);

    double s0 = ray->end()(0) - q0;
    double s1 = ray->end()(1) - q1;
    */
    Point bumped_incoming_point = incoming_point + 1e-9 * ray_direction; // + 1e-10 * ray_direction;

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
    {
      //      libMesh::err<<"bumped_incoming_point: "<<bumped_incoming_point<<std::endl;
      libMesh::err << "ray_direction: " << ray_direction << std::endl;
    }
#endif

    /*
    double q0 = incoming_point(0) + 1e-9*(ray->end()(0) - incoming_point(0));
    double q1 = incoming_point(1) + 1e-9*(ray->end()(1) - incoming_point(1));
    */

    double intersection_distance = 0;
    //    double centroid_distance = std::numeric_limits<Real>::max();

    Real t;

    // The x,y coordinates of this will be the bilinear coodinates of the intersection point
    Point u_v;

    // Center of the bilinear coordinates
    //    Point centroid(0.5, 0.5, 0);

    //    auto & normals = elem_normals.find(current_elem)->second;

    for (unsigned long int i = 0; i < n_sides; i++)
    {
      if (i == incoming_side) // Don't search backwards
        continue;

        // Backface culling
        //      if (normals[i] * ray_direction < -TOLERANCE)
        //        continue;

#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        libMesh::err << "Trying side: " << i << std::endl;
#endif

      bool intersected = /*intersectQuad*/
          intersectQuadUsingTriangles(bumped_incoming_point,
                                      ray_direction,
                                      *current_elem->get_node(Hex8::side_nodes_map[i][3]),
                                      *current_elem->get_node(Hex8::side_nodes_map[i][2]),
                                      *current_elem->get_node(Hex8::side_nodes_map[i][1]),
                                      *current_elem->get_node(Hex8::side_nodes_map[i][0]),
                                      u_v(0),
                                      u_v(1),
                                      t,
                                      ray);

      if (intersected)
      {
#ifdef USE_DEBUG_RAY
        if (DEBUG_IF)
          libMesh::err << "Intersected" << std::endl;
#endif
        if (t > intersection_distance)
        {
          //          Elem * neighbor = current_elem->neighbor(i);

#ifdef USE_DEBUG_RAY
          if (DEBUG_IF)
            libMesh::err << "Largest Intersection distance: " << t << std::endl;
#endif
          intersected_side = i;
          intersection_distance = t;

          intersection_point = incoming_point + t * ray_direction;
          /*
                    if (t > 1e-9)
                      return intersected_side;
          #ifdef USE_DEBUG_RAY
                    else if (DEBUG_IF) libMesh::err << t << " Not larger than 1e-9 so rejected " <<
          std::endl; #endif
          */
        }

        /*
        double current_centroid_distance = (centroid - u_v).norm();

        if (current_centroid_distance < centroid_distance)
        {
          Elem * neighbor = current_elem->neighbor(i);

          intersected_side = i;
          intersection_distance = t;
          centroid_distance = current_centroid_distance;

          intersection_point = incoming_point + t*ray_direction;
        }
        */
      }
    }
  }

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
  {
    libMesh::err << "Found side: " << intersected_side << std::endl;
    libMesh::err << "Intersection point: " << intersection_point << std::endl;
  }
#endif
  return intersected_side;
}

int
sideIntersectedByLine1D(const Elem * current_elem,
                        unsigned long int /*incoming_side*/,
                        const Point & /*incoming_point*/,
                        const std::shared_ptr<Ray> & /*ray*/,
                        const Point & in_ray_direction,
                        Point & intersection_point,
                        Point & /*boundaryintersection_point*/)
{
  Real ray_direction = in_ray_direction(0);

  // Ray is moving left: node/side 0 is the left node for EDGE elements
  if (ray_direction < 0)
  {
    intersection_point = *current_elem->get_node(0);
    return 0;
  }
  // Ray is moving right: node/side 1 is the right node for EDGE elements
  else
  {
    intersection_point = *current_elem->get_node(1);
    return 1;
  }

  return -1;
}

/**
 * Figure out which (if any) side of an Elem is intersected by a line.
 *
 * @param elem The elem to search
 * @param not_side Sides to _not_ search (Use -1 if you want to search all sides)
 * @param intersection_point If an intersection is found this will be filled with the x,y,z
 * position
 * of that intersection
 * @return The side that is intersected by the line.  Will return -1 if it doesn't intersect any
 * side
 */
/*
int sideIntersectedByLine(const Elem * elem, std::vector<int> & not_side, const Ray & ray, Point &
intersection_point)
{
  unsigned long int n_sides = elem->n_sides();

  // Whether or not they intersect
  bool intersect = false;

  unsigned long int dim = elem->dim();

  for (unsigned long int i=0; i<n_sides; i++)
  {
    // Don't search the "not_side"
    // Note: A linear search is fine here because this vector is going to be < n_sides
    if (std::find(not_side.begin(), not_side.end(), static_cast<int>(i)) != not_side.end())
      continue;

    // Get a simplified side element
    UniquePtr<Elem> side_elem = elem->side(i);

    if (dim == 3)
    {
      // Make a plane out of the first three nodes on the side
      Plane plane(side_elem->point(0), side_elem->point(1), side_elem->point(2));

      // See if they intersect
      intersect = ray.intersect(plane, intersection_point);
    }
    else if (dim == 2)
    {
      // Make a Line Segment out of the first two nodes on the side
      Ray side_segment(side_elem->point(0), side_elem->point(1));

      // See if they intersect
      intersect = ray.intersect(side_segment, intersection_point);
    }
    else // 1D
    {
      // See if the line segment contains the point
      intersect = ray.contains_point(side_elem->point(0));

      // If it does then save off that one point as the intersection point
      if (intersect)
        intersection_point = side_elem->point(0);
    }

    if (intersect)
    {
      if (side_elem->contains_point(intersection_point))
      {
        Elem * neighbor = elem->neighbor(i);

        // If this side is on a boundary, let's do another search and see if we can find a better
candidate
        if (!neighbor)
        {
          not_side.push_back(i); // Make sure we don't find this side again

          int better_side = sideIntersectedByLine(elem, not_side, ray, intersection_point);

          if (better_side != -1)
            return better_side;
        }

        return i;
      }
    }
  }

  // Didn't find one
  return -1;
}
*/


/**
 * Returns the side number for elem that neighbor is on
 *
 * Returns -1 if the neighbor can't be found to be a neighbor
 */
/*
int
sideNeighborIsOn(const Elem * elem, const Elem * neighbor)
{
 unsigned long int n_sides = elem->n_sides();

 for (unsigned long int i = 0; i < n_sides; i++)
 {
   if (elem->neighbor(i) == neighbor)
     return i;
 }

 return -1;
}
*/
void
TraceRay::find_point_neighbors(
    const Elem * current_elem,
    const Point & p,
    MooseUtils::StaticallyAllocatedSet<const Elem *, MAX_POINT_NEIGHBORS> & neighbor_set,
    const std::shared_ptr<Ray> &
#ifdef USE_DEBUG_RAY
        ray
#endif
)
{
  libmesh_assert(current_elem->contains_point(p));
  libmesh_assert(current_elem->active());

  neighbor_set.clear();
  neighbor_set.insert(current_elem);

  MooseUtils::StaticallyAllocatedSet<const Elem *, MAX_POINT_NEIGHBORS> untested_set;
  MooseUtils::StaticallyAllocatedSet<const Elem *, MAX_POINT_NEIGHBORS> next_untested_set;

  untested_set.insert(current_elem);

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "find_point_neighbors(" << current_elem->id() << ")" << std::endl;
#endif

  while (!untested_set.empty())
  {
    // Loop over all the elements in the patch that haven't already
    // been tested
    auto it = untested_set.begin();
    const auto end = untested_set.end();

    for (; it != end; ++it)
    {
      const Elem * elem = *it;

      const auto n_sides = elem->n_sides();

      for (unsigned long int s = 0; s < n_sides; s++)
      {
        const Elem * current_neighbor = elem->neighbor_ptr(s);
#ifdef USE_DEBUG_RAY
        if (DEBUG_IF && current_neighbor)
          std::cerr << " Testing neighbor " << current_neighbor->id() << std::endl;

        if (DEBUG_IF && current_neighbor == remote_elem)
          std::cerr << " Neighbor is remote! " << std::endl;
#endif

        if (current_neighbor &&
            current_neighbor != remote_elem) // we have a real neighbor on current_elem side
        {
          if (current_neighbor->active()) // ... if it is active
          {
            if (current_neighbor->close_to_point(p,
                                                 1e-5) /*contains_point(p)*/) // ... and touches p
            {
              // Make sure we'll test it
              if (!neighbor_set.contains(current_neighbor))
                next_untested_set.insert(current_neighbor);

#ifdef USE_DEBUG_RAY
              if (DEBUG_IF)
                std::cerr << "  Contained in " << current_neighbor->id() << std::endl;
#endif
              // And add it
              neighbor_set.insert(current_neighbor);
            }
          }
#ifdef LIBMESH_ENABLE_AMR
          else // ... the neighbor is *not* active,
          {    // ... so add *all* neighboring
               // active children that touch p
            current_neighbor->active_family_tree_by_neighbor(_active_neighbor_children, elem);

            std::vector<const Elem *>::const_iterator child_it = _active_neighbor_children.begin();
            const std::vector<const Elem *>::const_iterator child_end =
                _active_neighbor_children.end();
            for (; child_it != child_end; ++child_it)
            {
              const Elem * current_child = *child_it;
              if (current_child->contains_point(p))
              {
                // Make sure we'll test it
                if (!neighbor_set.contains(current_child))
                  next_untested_set.insert(current_child);

                neighbor_set.insert(current_child);
              }
            }
          }
#endif // #ifdef LIBMESH_ENABLE_AMR
        }
      }
    }
    untested_set.swap(next_untested_set);
    next_untested_set.clear();
  }
}

Elem *
TraceRay::recursivelyFindChildContainingPoint(const Elem * current_child,
                                              const Point & intersection_point,
                                              const std::vector<unsigned long int> & children_ids)
{
  for (auto child_id : children_ids)
  {
    auto child = current_child->child(child_id);

    if (child->contains_point(intersection_point))
    {
      if (child->active())
        return child;
      else
        return recursivelyFindChildContainingPoint(child, intersection_point, children_ids);
    }
  }

  std::ostringstream children;

  for (auto child_id : children_ids)
  {
    auto child = current_child->child(child_id);

    children << child_id << ":\n" << *child << std::endl;
  }

  mooseError("Unable to find child containing the point!  Should never happen! Elem:",
             *current_child,
             " intersection_point: ",
             intersection_point,
             "\nchildren searched: \n",
             children.str());
}

Elem *
TraceRay::childOnSide(const Elem * current_elem, const Point & p, unsigned long int side)
{
  // Get the children on that side
  auto & children_ids = quad4_side_to_children[side];

  return recursivelyFindChildContainingPoint(current_elem, p, children_ids);
}

Elem *
TraceRay::getNeighbor(const Elem * current_elem,
                      unsigned long int intersected_side,
                      Point & intersection_point)
{
  return current_elem->neighbor(intersected_side);

  auto neighbor = current_elem->neighbor(intersected_side);

  if (!neighbor || neighbor->active())
    return neighbor;

  else // There is adaptivity... need to find the active child that contains the point
  {
    // Get the side the current elem occupies for the neighbor
    auto neighbor_side =
        neighbor->which_neighbor_am_i(current_elem); // sideNeighborIsOn(neighbor, current_elem);

    return childOnSide(neighbor, intersection_point, neighbor_side);
  }
}

/**
 * Checks to see if the ray is perfectly hitting any nodes on this element
 */
void
TraceRay::possiblyHittingNode(const std::shared_ptr<Ray> &

#ifdef USE_DEBUG_RAY
                                  ray
#endif
                              ,
                              const Point & incoming_point,
                              const Elem * current_elem,
                              const Point & ray_direction,
                              int & node_hit,
                              Point & intersection_point,
                              int & intersected_side)
{
  // To do this we're simply going to see if the direction from the incoming point
  // to each node is the same as the direction the Ray is traveling in... simple

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "Checking to see if the ray is perfectly hitting a node:" << std::endl;
#endif

  for (unsigned int n = 0; n < current_elem->n_nodes(); n++)
  {
    const auto & node_ptr = current_elem->node_ptr(n);

    // Make a normalized direction vector
    Point node_direction = (*node_ptr) - incoming_point;

    auto distance = node_direction.size();

    // Don't collide with a node we're already at
    if (distance < 1e-8)
      continue;

    node_direction /= distance;

    // Now see if things are in the same direction

    //    std::cout << "Node direction: " << node_direction << std::endl;
    //    std::cout << "Ray direction: " << ray_direction << std::endl;

    auto dot = node_direction * ray_direction;

    //    if (dot > 0.8)
    //      std::cout << "Dot: " << dot << std::endl;
#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
      std::cerr << " " << n << " Node Dot: " << dot << std::endl;
#endif

    // Did we find a match?
    // This number was experimentally determined by slowly ramping it up
    // Might need further tweaking
    if (dot > 0.99999999)
    {
      node_hit = n;

      std::cerr << "Hit node! Dot: " << dot << std::endl;

      // We're going to intersect at the node
      intersection_point = *node_ptr;

      // Pick any side the node is on
      for (unsigned int s = 0; s < current_elem->n_sides(); s++)
      {
        if (current_elem->is_node_on_side(n, s))
        {
          intersected_side = s;
          break;
        }
      }

      break;
    }
  }
}

/**
 * Checks to see if the ray will hit a "edge" (an edge of an element in 3D)
 */
void
TraceRay::possiblyHittingEdge(const std::shared_ptr<Ray> & ray,
                                const Point & incoming_point,
                                const Elem * current_elem,
                                const Point & ray_direction,
                                int & edge_hit,
                                Point & intersection_point,
                                int & intersected_side)
{
  if (current_elem->type() != HEX8)
  {
    std::cerr << "possiblyHittingEdge() is only implemented for HEX8 elements!" << std::endl;
    std::abort();
  }

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
      std::cerr << "Checking element edges: " << std::endl;
#endif
  Real intersection_distance = 0;

  Real best_intersection_distance = 1e-8;
  Point best_intersection_point;
  int best_side = -1;
  int best_edge = -1;

  // Loop over all of the edges and test them
  for (unsigned int i = 0; i < 12; i++)
  {
    auto & edge_nodes = Hex8::edge_nodes_map[i];

    auto & n0 = current_elem->node_ref(edge_nodes[0]);
    auto & n1 = current_elem->node_ref(edge_nodes[1]);

    bool intersected = lineLineIntersect3DVanilla(ray, incoming_point, ray->end(), n0, n1, intersection_point, intersection_distance);

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
      std::cerr << " Edge" << i << ": " << intersected << std::endl;
#endif

    if (intersected && intersection_distance > best_intersection_distance)
    {
#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
      {
        std::cerr << "  best_intersection_distance: " << intersection_distance << std::endl;
        std::cerr << "  best_intersection_point: " << intersection_point << std::endl;
        std::cerr << "  best_edge: " << i << std::endl;
        std::cerr << "  best_side: " << hex8_edge_to_sides[i][0] << std::endl;
      }
#endif

      best_intersection_distance = intersection_distance;
      best_intersection_point = intersection_point;
      best_edge = i;
      best_side = hex8_edge_to_sides[i][0];
    }
  }

  if (best_edge != -1)
  {
    intersection_point = best_intersection_point;
    edge_hit = best_edge;
    intersected_side = best_side;
  }
}


void
TraceRay::endPossiblyOnBoundarySide(const std::shared_ptr<Ray> & ray,
                                    const Elem * current_elem,
                                    Point & intersection_point,
                                    int & intersected_side)
{
  const auto n_sides = current_elem->n_sides();

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "Does the endpoint lie within a boundary side?" << std::endl;
#endif

  // Loop over the sides and see if they are on the boundary - if they are, then
  // see if they contain the end point
  for (auto s = 0u; s < n_sides; s++)
  {
    if (current_elem->neighbor(s) == NULL)
    {
      auto side_elem = current_elem->build_side(s);

#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        std::cerr << " Checking boundary side: " << s << std::endl;
#endif

      if (side_elem->contains_point(ray->end()))
      {
#ifdef USE_DEBUG_RAY
        if (DEBUG_IF)
          std::cerr << "  Contained in side!" << std::endl;
#endif

        intersection_point = ray->end();
        intersected_side = s;
        break;
      }
    }
  }
}

/**
 * This routine checks to see if the current start point of the ray is already REALLY
 * close to the side of the domeain where it will need to be stopped.  This keeps badness
 * from happening when a ray strikes a side right at a junction that is on the boundary -
 * but happens to pick a side that is not _actually_ on the boundary.
 */
void
TraceRay::possiblyOnBoundary(const std::shared_ptr<Ray> & ray,
                             const Point & incoming_point,
                             const Elem * current_elem,
                             unsigned long int incoming_side,
                             Point & intersection_point,
                             int & intersected_side)
{
  if ((ray->end() - incoming_point).size() < 1e-8) // Only allow this if we're near the end
  {
    const auto n_sides = current_elem->n_sides();

    unsigned long int sides_on_boundary = 0;
    // First, see how many sides are by the boundary
    for (auto s = 0u; s < n_sides; s++)
      if (current_elem->neighbor(s) == NULL)
        sides_on_boundary++;

    // Can't find a boundary side if this element has none
    if (!sides_on_boundary)
      return;

    // See if we're on the edge of the domain
    _work_point = incoming_point - _b_box.min();
    _work_point2 = incoming_point - _b_box.max();

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
    {
      std::cerr << "Work point: " << _work_point << std::endl;
      std::cerr << "Work point2: " << _work_point2 << std::endl;
    }
#endif

    // We're going to look at our bounding box for the domain to know if we're
    // close to the side.  What we're looking for here are dims where the difference
    // between the Ray's dim and the bounding box dim is zero
    unsigned long int num_zero = 0;

    for (unsigned long int i = 0; i < _mesh_dim; i++)
      if (MooseUtils::absoluteFuzzyEqual(std::abs(_work_point(i)), 0., 5e-5))
        num_zero++;

    for (unsigned long int i = 0; i < _mesh_dim; i++)
      if (MooseUtils::absoluteFuzzyEqual(std::abs(_work_point2(i)), 0., 5e-5))
        num_zero++;

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
      std::cerr << "Checking to see if by boundary: " << num_zero << std::endl;
#endif

    if (num_zero) // If we are, then let's just call this point good
    {
#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        std::cerr << "On domain boundary! And have moved: " << ray->distance() << std::endl;
#endif

      // Now what we're going to do is search the sides of this element for a side
      // on the boundary that still contains this point
      if (sides_on_boundary == 1) // If there's only one let's just pick it
      {
        for (auto s = 0u; s < n_sides; s++)
        {
          if (current_elem->neighbor(s) == NULL)
          {
            intersected_side = s;
            intersection_point = ray->end();

            return;
          }
        }
      }
      else if (sides_on_boundary > 1) // If there's more than one then we need to choose
      {
        for (auto s = 0u; s < n_sides; s++)
        {
          if (s != incoming_side) // Don't allow a ray to turn around!
          {
            // Only check sides that are up against the domain edge
            if (current_elem->neighbor(s) == NULL)
            {
#ifdef USE_DEBUG_RAY
              if (DEBUG_IF)
                std::cerr << "Checking to see if side " << s << "contains the point" << std::endl;
#endif
              auto side_elem = current_elem->build_side(s);

              if (side_elem->contains_point(incoming_point))
              {
#ifdef USE_DEBUG_RAY
                if (DEBUG_IF)
                {
                  std::cerr << ray_count << " Found boundary side: " << s << std::endl;
                  std::cerr << ray_count << " Boundary intersection point: " << incoming_point
                            << std::endl;
                }
#endif
                intersected_side = s;
                intersection_point = incoming_point;
                break;
              }
            }
          }
        }
      }
    }
  }
}

void
TraceRay::tryToMoveThroughPointNeighbors(const Elem * current_elem,
                                         const ElemType elem_type,
                                         const Point & incoming_point,
                                         const std::shared_ptr<Ray> & ray,
                                         const Point & ray_direction,
                                         Point & boundary_intersection_point,
                                         int & intersected_side,
                                         const Elem *& best_neighbor,
                                         unsigned long int & best_side)
{
  // Let's first try grabbing the point_neighbors for this element to see if they are good
  // candidates
  find_point_neighbors(current_elem, incoming_point, _point_neighbors, ray);

  Point intersection_point;

#ifdef USE_DEBUG_RAY
  if (DEBUG_IF)
    std::cerr << "Num point_neighbors: " << _point_neighbors.size() << std::endl;
#endif

  Real longest_distance = 0;

  // Try the search on each neighbor and see if something good happens
  for (auto & neighbor : _point_neighbors)
  {
    if (neighbor != current_elem) // Don't need to look through this element again
    {
#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        std::cerr << "Trying neighbor " << neighbor->id() << std::endl;
#endif

      // First find the side the point is on in the neighbor
      auto side = -1;

      const auto n_sides = neighbor->n_sides();

      for (auto s = 0u; s < n_sides; s++)
      {
        auto side_elem = neighbor->build_side(s);

        if (side_elem->contains_point(incoming_point))
        {
          side = s;
          break;
        }
      }

      if (side == -1)
      {
#ifdef USE_DEBUG_RAY
        if (DEBUG_IF)
          std::cerr << "Neighbor sides didn't contain the point" << neighbor->id() << std::endl;
#endif
        continue;
      }

      // Now try to find an intersection in that neighbor
      if (elem_type == QUAD4)
        intersected_side = sideIntersectedByLine2D<Quad4>(neighbor,
                                                          side,
                                                          incoming_point,
                                                          ray,
                                                          ray_direction,
                                                          _elem_normals,
                                                          intersection_point,
                                                          boundary_intersection_point);
      else if (elem_type == TRI3)
        intersected_side = sideIntersectedByLine2D<Tri3>(neighbor,
                                                         side,
                                                         incoming_point,
                                                         ray,
                                                         ray_direction,
                                                         _elem_normals,
                                                         intersection_point,
                                                         boundary_intersection_point);
      else if (elem_type == HEX8)
        intersected_side = sideIntersectedByLineHex8(neighbor,
                                                     side,
                                                     incoming_point,
                                                     ray,
                                                     ray_direction,
                                                     _elem_normals,
                                                     intersection_point,
                                                     boundary_intersection_point);

      // If we found an intersection let's go with it
      if (intersected_side != -1)
      {
        auto distance = (intersection_point - incoming_point).norm();

        if (distance > longest_distance)
        {
          best_neighbor = neighbor;
          best_side = side;
          longest_distance = distance;
        }
      }
    }
  }
}

void
TraceRay::checkForCornerHitAndApplyBCs(const Elem * current_elem,
                                       const Point & intersection_point,
                                       std::shared_ptr<Ray> & ray)
{
  // Need to see if the ray hit _right_ on the "corner" of the domain
  // If it did, then we're going to apply all of the boundary conditions for each
  // side at that corner

  // First determine if we hit a corner.
  // Subtract off the current point from the min/max of the bounding box
  // If any two entries from that subtraction are (near) zero then we're at a corner
  _work_point = intersection_point - _b_box.min();
  _work_point2 = intersection_point - _b_box.max();

  unsigned long int num_zero = 0;

  for (unsigned long int i = 0; i < _mesh_dim; i++)
    if (MooseUtils::absoluteFuzzyEqual(std::abs(_work_point(i)), 0., 1e-6))
      num_zero++;

  for (unsigned long int i = 0; i < _mesh_dim; i++)
    if (MooseUtils::absoluteFuzzyEqual(std::abs(_work_point2(i)), 0., 1e-6))
      num_zero++;

  //        if (ray->id() == 4811)
  //          libMesh::err<<"Here!"<<std::endl;

  if (num_zero > 1)
  {
    //          if (ray->id() == 4811)
    //          libMesh::err<<ray_count<<" Ray hit a domain corner!
    //          "<<intersection_point<<std::endl;

    find_point_neighbors(current_elem, intersection_point, _point_neighbors, ray);

    // Try the search on each neighbor and see if something good happens
    for (auto & neighbor : _point_neighbors)
    {
//            if (neighbor != current_elem) // Don't need to look through this element again
//            {
#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        libMesh::err << "Looking at neighbor: " << neighbor->id() << std::endl;
#endif

      // First find the side the point is on in the neighbor
      auto side = -1;

      const auto n_sides = neighbor->n_sides();

      for (auto s = 0u; s < n_sides; s++)
      {
        auto side_elem = neighbor->build_side(s);

        if (side_elem->contains_point(intersection_point))
        {
          side = s;

          _mesh.get_boundary_info().boundary_ids(neighbor, side, _ids);

          for (const auto & bid : _ids)
          {
            // Don't apply the same BC twice!
            if (!_applied_ids.contains(bid))
            {
              _applied_ids.insert(bid);
              onBoundary(current_elem, intersection_point, side, bid, neighbor, ray);
            }
          }
        }
      }
    }
  }
}

/**
 */
void
TraceRay::trace(std::shared_ptr<Ray> & ray)
{
  //  TIME_SECTION(_trace_timer);

  const Elem * current_elem = NULL;
  int intersected_side = -1;
  Point intersection_point;
  Point boundary_intersection_point;

  Point ray_direction = (ray->end() - ray->start());
  ray_direction /= ray_direction.norm();

  Point incoming_point = ray->start();

  current_elem = ray->startingElem();

  unsigned long int incoming_side = ray->incomingSide();

  // The node number of the current element that might have been hit by the ray
  // -1 means it didn't hit a node on the current element
  int node_hit = -1;

  // The edge number of the current element that might have been hit by the ray
  // -1 means it didn't hit a edge on the current element
  int edge_hit = -1;

  bool keep_tracing = true;

  unsigned long int count = 0;

  unsigned long int halo_count = 0;

  ray_count++;

  SubdomainID current_subdomain = -1;

#ifdef USE_DEBUG_RAY
  Mesh * debug_mesh = NULL;

  unsigned long int debug_node_count = 0;
#endif

  auto pid = _mesh.comm().rank();

  //#ifdef USE_DEBUG_RAY
  //  if (pid != debug_ray_pid)
  //    sleep(1000);
  //#endif

  _mesh_dim = _mesh.mesh_dimension();

#ifdef USE_DEBUG_RAY
  if (false && DEBUG_IF && pid == 0)
  //  if (pid == 1)
  {
    debug_mesh = new Mesh(Parallel::Communicator(), _mesh.mesh_dimension());
    debug_mesh->skip_partitioning(true);
  }
#endif

//  libMesh::err<<std::endl<<"Starting new ray!"<<std::endl;
#ifdef USE_DEBUG_RAY
  // Add an element for the ray

  if (debug_mesh && DEBUG_IF && pid == 0)
  //  if (pid == 1)
  {

    debug_mesh->add_point(ray->start());
    debug_node_count++;
    debug_mesh->add_point(ray->end());
    debug_node_count++;

    Elem * elem = Elem::build(EDGE2).release();
    elem->subdomain_id() = 1;

    elem = debug_mesh->add_elem(elem);

    elem->set_node(0) = debug_mesh->node_ptr((debug_node_count - 2));
    elem->set_node(1) = debug_mesh->node_ptr((debug_node_count - 2) + 1);
  }

#endif

  do
  { // Find the side of this element that the Ray intersects... while ignoring the incoming side
    // (we
    // don't want to move backward!)
    intersected_side = -1;

    node_hit = -1;
    edge_hit = -1;

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
    {
      std::cerr << "\n\n" << pid << " At top!" << std::endl;
      std::cerr << pid << " Debug ray UID: " << ray->id() << std::endl;
      std::cerr << pid << " Ray Start: " << ray->start() << std::endl;
      std::cerr << pid << " Ray End: " << ray->end() << std::endl;
      std::cerr << pid << " Traveled: " << ray->distance() << std::endl;
      std::cerr << pid << " To Go: " << (ray->end() - ray->start()).size() << std::endl;
    }
#endif

#ifdef USE_DEBUG_RAY
    if (DEBUG_IF)
    {
      std::cerr << pid << " current_elem: " << current_elem->id() << std::endl;
      std::cerr << pid << *current_elem << std::endl;
    }
#endif

    bool ends_in_elem = false;

    // If the ray ends within the mesh... let's see if it ends within this element!
    if (ray->endsWithinMesh())
    {
#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        std::cerr << pid << " ends within mesh!" << std::endl;
#endif
      // Check the ending elem id
      if (ray->endingElemId() == current_elem->id())
        ends_in_elem = true;
      // Ending elem id is not set: check with contains_point
      else if (ray->endingElemId() == DofObject::invalid_id &&
               current_elem->contains_point(ray->end()))
        ends_in_elem = true;

      if (ends_in_elem)
      {
#ifdef USE_DEBUG_RAY
        if (DEBUG_IF)
          std::cerr << pid << " ends within current elem!" << std::endl;
#endif
        intersection_point = ray->end();
      }
    }

    auto elem_type = current_elem->type();

    if (!ends_in_elem)
    {
      //      TIME_SECTION(_intersection_timer);

      /*

        libMesh::err<<"In element: "<<current_elem->id()<<std::endl;

    if (debug_mesh && ray->id() == debug_ray_id && pid == debug_ray_pid)
  //  if ( == 1)
      {
        for (unsigned long int i=0; i<current_elem->n_nodes(); i++)
        {
          debug_mesh->add_point(*current_elem->get_node(i));
          debug_node_count++;
        }

  //      Elem * debug_elem = Elem::build(QUAD4).release();
        Elem * debug_elem = Elem::build(HEX8).release();
        debug_elem->subdomain_id() = 2;

        debug_elem = debug_mesh->add_elem(debug_elem);

        for (unsigned long int i=0; i<current_elem->n_nodes(); i++)
        {
  //        libMesh::err<<"Elem nodes: "<<*debug_mesh->node_ptr((debug_node_count -
  current_elem->n_nodes()) + i)<<std::endl;
          debug_elem->set_node(i) = debug_mesh->node_ptr((debug_node_count -
  current_elem->n_nodes())
  + i);
        }

      }


    if (ray->id() == debug_ray_id && pid == debug_ray_pid)
    {
      for (unsigned long int s=0; s<current_elem->n_sides(); s++)
      {
        fe_face->reinit(current_elem, s);
        libMesh::err<<"Side "<<s<<" normal: "<<normals[0]<<std::endl;
      }
    }

      */
      if (elem_type == QUAD4)
        intersected_side = sideIntersectedByLine2D<Quad4>(current_elem,
                                                          incoming_side,
                                                          incoming_point,
                                                          ray,
                                                          ray_direction,
                                                          _elem_normals,
                                                          intersection_point,
                                                          boundary_intersection_point);
      else if (elem_type == TRI3)
        intersected_side = sideIntersectedByLine2D<Tri3>(current_elem,
                                                         incoming_side,
                                                         incoming_point,
                                                         ray,
                                                         ray_direction,
                                                         _elem_normals,
                                                         intersection_point,
                                                         boundary_intersection_point);
      else if (elem_type == HEX8)
        intersected_side = sideIntersectedByLineHex8(current_elem,
                                                     incoming_side,
                                                     incoming_point,
                                                     ray,
                                                     ray_direction,
                                                     _elem_normals,
                                                     intersection_point,
                                                     boundary_intersection_point);

      else if (elem_type == EDGE2 || elem_type == EDGE3 || elem_type == EDGE4)
        intersected_side = sideIntersectedByLine1D(current_elem,
                                                   incoming_side,
                                                   incoming_point,
                                                   ray,
                                                   ray_direction,
                                                   intersection_point,
                                                   boundary_intersection_point);
    }

    if (intersected_side == -1 && !ends_in_elem) // If we failed to find a side... try harder
    {
      //      TIME_SECTION(_try_harder_timer);

#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
      {
        std::cerr << pid << " Didn't find a side... trying harder" << std::endl;
        std::cerr << pid << " ray->distance(): " << ray->distance() << std::endl;
      }
#endif

      // Check to see if the ray perfectly strikes any nodes of this element
      if (intersected_side == -1)
        possiblyHittingNode(ray,
                            incoming_point,
                            current_elem,
                            ray_direction,
                            node_hit,
                            intersection_point,
                            intersected_side);

      if (intersected_side == -1 && _mesh_dim == 3)
        possiblyHittingEdge(ray,
                            incoming_point,
                            current_elem,
                            ray_direction,
                            edge_hit,
                            intersection_point,
                            intersected_side);

      if (intersected_side == -1)
        endPossiblyOnBoundarySide(ray, current_elem, intersection_point, intersected_side);

      if (intersected_side == -1)
        possiblyOnBoundary(
            ray, incoming_point, current_elem, incoming_side, intersection_point, intersected_side);

      if (intersected_side == -1) // Need to do a more exhaustive search
      {
        const Elem * best_neighbor = NULL;
        unsigned long int best_side = 0;

        tryToMoveThroughPointNeighbors(current_elem,
                                       elem_type,
                                       incoming_point,
                                       ray,
                                       ray_direction,
                                       boundary_intersection_point,
                                       intersected_side,
                                       best_neighbor,
                                       best_side);

        if (best_neighbor)
        {
#ifdef USE_DEBUG_RAY
          if (DEBUG_IF)
            std::cerr << pid << " best_neighbor: " << best_neighbor->id()
                      << " best_side: " << best_side << std::endl;
#endif
          current_elem = best_neighbor;
          incoming_side = best_side;

          if (current_elem->processor_id() != pid)
          {
#ifdef USE_DEBUG_RAY
            if (DEBUG_IF)
              std::cerr << pid << " point neighbor off processor: " << current_elem->processor_id()
                        << std::endl;
#endif
            ray->setStartingElem(current_elem);
            ray->setIncomingSide(incoming_side);
            ray->setStart(incoming_point);
            break;
          }

          continue;
        }
        else
        {

#ifdef USE_DEBUG_RAY
          if (DEBUG_IF)
#endif
          if (!_tolerate_failure)
          {
            std::cerr << pid << " Ray " << ray->id()
                      << " is unable to find a point neighbor to move through!" << std::endl;

            std::cerr << pid << " All options have been exhausted... failing!" << std::endl;

            std::abort();
          }
        }



        count++;

        if (count > 200000)
        {

#ifdef USE_DEBUG_RAY
          if (debug_mesh && DEBUG_IF)
          //  if (pid == 1)
          {

            debug_mesh->prepare_for_use();
            debug_mesh->write("debug.e");
          }
#endif

#ifdef USE_DEBUG_RAY
          if (DEBUG_IF)
#endif
          libMesh::err << pid << " "
                       << "ray_count: " << ray_count << " ray_id: " << ray->id()
                       << " Endless loop while trying harder!\n"
                       << " Start: " << ray->start() << "\n"
                       << " End: " << ray->end() << " Distance: " << ray->distance() << std::endl;

          ray->setShouldContinue(false);

#ifdef USE_DEBUG_RAY
          if (DEBUG_IF)
#endif
          if (!_tolerate_failure)
            std::abort();

          break;
        }

        // This causes the loop to restart with the neighbor as current_elem
        continue;
      }

      /*
      if (intersected_side == -1)
      {
        std::vector<int> not_side = {incoming_side};

        intersected_side = sideIntersectedByLine(current_elem, not_side, ray, intersection_point);
      }
      */
    }

    /*
    #ifdef USE_DEBUG_RAY
      if (debug_mesh && ray->id() == debug_ray_id && pid == debug_ray_pid)
      {
        libMesh::err<<"Adding segment:"<<std::endl;
        libMesh::err<<"  "<<incoming_point<<std::endl;
        libMesh::err<<"  "<<intersection_point<<std::endl;

        debug_mesh->add_point(incoming_point);
        debug_node_count++;
        debug_mesh->add_point(intersection_point);
        debug_node_count++;

        Elem * elem = Elem::build(EDGE2).release();
        elem->subdomain_id() = 1;

        elem = debug_mesh->add_elem(elem);

        elem->set_node(0) = debug_mesh->node_ptr((debug_node_count - 2));
        elem->set_node(1) = debug_mesh->node_ptr((debug_node_count - 2) + 1);

    //    debug_mesh->prepare_for_use();
    //    debug_mesh->write("debug.e");

        }
    #endif
    */

    if (intersected_side != -1 || ends_in_elem) // This is true if we have found an intersection
    {
      //      TIME_SECTION(_found_intersection);

#ifdef USE_DEBUG_RAY
      if (DEBUG_IF)
        std::cerr << pid << " Found a side!" << std::endl;
#endif

      ray->intersections()++;
      _num_intersections++;

      //      libMesh::err<<"Intersected side "<<intersected_side<<" at
      //      "<<intersection_point<<std::endl;

      // Use a work point here so temporaries don't get created
      _work_point = intersection_point;
      _work_point -= incoming_point;

      // The 0.79... is for TY quadrature in 2D
      // TODO: fix this better for 3D
      ray->distance() += _work_point.norm();
      /*
            // Need to reinit at a physical point.  For now, we'll just take that to be halfway
         along the line
            // TODO: In the future we'll use quadrature here to integrate along the line
            phys_points[0] = (intersection_point + incoming_point) / 2.0;

            fe_problem.reinitElemPhys(current_elem, phys_points, tid);
      */

      if (current_elem->subdomain_id() != current_subdomain)
      {
        current_subdomain = current_elem->subdomain_id();

        subdomainSetup(current_subdomain, ray);
      }

      onSegment(current_elem, incoming_point, intersection_point, ends_in_elem);

      // Check if a kernel set to kill the ray
      if (!ray->shouldContinue())
        break;

      // If the Ray is beyond the max distance we need to kill it
      if (ray->distance() >= _ray_max_distance)
      {
#ifdef USE_DEBUG_RAY
        if (DEBUG_IF)
          std::cerr << pid << " Killing because at max distance!" << std::endl;
#endif
        //        libMesh::err<<"Ray over max distance: "<<ray->distance()<<std::endl;
        ray->setShouldContinue(false);
        break;
      }

      if (ends_in_elem)
      {
        ray->setShouldContinue(false);
        break;
      }

      // Grab the neighbor and set up where to go next
      {
        const Elem * neighbor = NULL;

        if (node_hit != -1 || edge_hit != -1)
        {
          // Since we hit a node - we know we're going to have trouble moving through the next
          // element
          // Let's go ahead and search for the best element to move through next

          // First, see if this node is on the boundary - if it is, then let's just use a null
          // neighbor
          bool on_boundary = false;

          if (node_hit != -1)
          {
            for (unsigned int s = 0; s < current_elem->n_sides(); s++)
            {
              if (current_elem->is_node_on_side(node_hit, s))
              {
                if (!current_elem->neighbor(s))
                {
                  on_boundary = true;
                  intersected_side = s;
//                  intersection_point = ray->end();
                  break;
                }
              }
            }
          }

          if (edge_hit != -1)
          {
            auto & sides = hex8_edge_to_sides[edge_hit];

            for (auto & side : sides)
            {
              if (!current_elem->neighbor(side))
              {
                on_boundary = true;
                intersected_side = side;
//                intersection_point = ray->end();
              }
            }
          }

          if (!on_boundary)
          {
            const Elem * best_neighbor = NULL;
            unsigned long int best_side = 0;

            tryToMoveThroughPointNeighbors(current_elem,
                                           elem_type,
                                           intersection_point,
                                           ray,
                                           ray_direction,
                                           boundary_intersection_point,
                                           intersected_side,
                                           best_neighbor,
                                           best_side);

            if (best_neighbor)
            {
#ifdef USE_DEBUG_RAY
              if (DEBUG_IF)
                std::cerr << pid << " best_neighbor: " << best_neighbor->id()
                          << " best_side: " << best_side << std::endl;
#endif
              current_elem = neighbor = best_neighbor;
              incoming_side = best_side;
              incoming_point = intersection_point;
            }
          }
        }
        else
        {
          // Get the neighbor on that side
          neighbor = getNeighbor(current_elem, intersected_side, intersection_point);

#ifdef USE_DEBUG_RAY
          if (DEBUG_IF)
            std::cerr << pid << " Neighbor: " << neighbor << std::endl;
#endif

          if (neighbor)
          {
#ifdef USE_DEBUG_RAY
            if (DEBUG_IF)
              std::cerr << pid << " Neighbor_id: " << neighbor->id() << std::endl;
#endif
            // Note: This is finding the side the current_elem is on for the neighbor.  That's the
            // "incoming_side" for the neighbor
            incoming_side = neighbor->which_neighbor_am_i(
                current_elem); // sideNeighborIsOn(neighbor, current_elem);

            // Recurse
            current_elem = neighbor;
          }
        }

        incoming_point = intersection_point;
        ray->setStart(incoming_point);

        // If the neighbor is off-processor... we need to set up the Ray data to be passed to
        // the
        // next processor
        if (neighbor && neighbor->processor_id() != pid)
        {
          halo_count++;

          if (halo_count == _halo_size)
          {
#ifdef USE_DEBUG_RAY
            if (DEBUG_IF)
              std::cerr << pid << " Ray going off processor!" << std::endl;
#endif
            ray->setStartingElem(neighbor);
            ray->setIncomingSide(incoming_side);
            ray->setStart(incoming_point);
            break;
          }
        }

        if (!neighbor) // Edge of the domain
        {
#ifdef USE_DEBUG_RAY
          if (DEBUG_IF)
            std::cerr << "At edge of domain!" << std::endl;
#endif

          /*

                    libMesh::err<<"Reflecting!"<<std::endl;
          */

          //        libMesh::err<<"Original:\n Start: "<<ray->start()<<"\n End:
          //        "<<ray->end()<<std::endl;

          // Need to see if this is a reflective BC
          // Get the "direction" of the Ray
          // Use _work_point so no temporaries are created

          //        libMesh::err<<" Direction: "<<_work_point<<std::endl;

          _mesh.get_boundary_info().boundary_ids(current_elem, intersected_side, _ids);

          for (const auto & bid : _ids)
          {
            if (!_applied_ids.contains(bid))
            {
              _applied_ids.insert(bid);
              onBoundary(current_elem, intersection_point, intersected_side, bid, NULL, ray);
            }
          }

          checkForCornerHitAndApplyBCs(current_elem, intersection_point, ray);

          finishedBoundary();

          _applied_ids.clear();

          // See if Ray was killed by the BC
          if (!ray->shouldContinue())
            break;
          else // TODO: This isn't quite right - should really read what was set by the BC
               // later...
          // but this will work for reflective for now
          {
            incoming_point = intersection_point;
            ray->setStart(incoming_point);

            incoming_side = intersected_side;
          }

          //        libMesh::err<<"Normal: "<<normals[0]<<std::endl;

          /*
      {

      debug_mesh->add_point(ray->start());
      debug_node_count++;
      debug_mesh->add_point(ray->end());
      debug_node_count++;

      Elem * elem = Elem::build(EDGE2).release();
      elem->subdomain_id() = 1;

      elem = debug_mesh->add_elem(elem);

      elem->set_node(0) = debug_mesh->node_ptr((debug_node_count - 2));
      elem->set_node(1) = debug_mesh->node_ptr((debug_node_count - 2) + 1);
      }

          */

          //        libMesh::err<<"Reflected:\n Start: "<<ray->start()<<"\n End:
          //        "<<ray->end()<<std::endl;
          //        libMesh::err<<" Direction: "<<_work_point<<std::endl;

          // ray->setStartingElem(NULL);
          // keep_tracing = false;
        }
      }
    }
    else // Shouldn't ever happen... but we'll just kill this Ray off.  It won't effect anything
    {

      /*
      #ifdef USE_DEBUG_RAY
        if (debug_mesh && ray->id() == debug_ray_id && pid == debug_ray_pid)
      //  if (pid == 1)
        {



          debug_mesh->add_point(incoming_point);
          debug_node_count++;
          debug_mesh->add_point(ray->end());
          debug_node_count++;

          Elem * elem = Elem::build(EDGE2).release();
          elem->subdomain_id() = 1;

          elem = debug_mesh->add_elem(elem);

          elem->set_node(0) = debug_mesh->node_ptr((debug_node_count - 2));
          elem->set_node(1) = debug_mesh->node_ptr((debug_node_count - 2) + 1);




          debug_mesh->prepare_for_use();
          debug_mesh->write("debug.e");
        }
      #endif
      */

      //      mooseError(ray_count << " Shouldn't be here\n" << ray->start() << "\n" <<
      //      ray->end());

      std::cerr << "Ray " << ray_count << " killed prematurely! "
                << " Unique ID: " << ray->id() << std::endl;
      std::cerr << "STUUUUFFFF" << std::endl;

#ifdef USE_DEBUG_RAY
      if (debug_mesh && DEBUG_IF)
      //  if (pid == 1)
      {

        debug_mesh->prepare_for_use();
        debug_mesh->write("debug.e");
      }
#endif

      if (!_tolerate_failure)
        std::abort();

      ray->setShouldContinue(false);
      keep_tracing = false;
    }

    count++;

    if (count > 2000000)
    {

#ifdef USE_DEBUG_RAY
      if (debug_mesh && DEBUG_IF)
      //  if (pid == 1)
      {

        debug_mesh->prepare_for_use();
        debug_mesh->write("debug.e");
      }
#endif
      libMesh::err << pid << " "
                   << "ray_count: " << ray_count << " ray_id: " << ray->id() << " Endless loop!"
                   << std::endl;
      libMesh::err << "STUUUUFFFF" << std::endl;
      ray->setShouldContinue(false);

      if (!_tolerate_failure)
        std::abort();

      break;
    }
  } while (keep_tracing);
}
