// Local Includes
#include "RayKernel.h"
#include "RayBC.h"
#include "RayProblem.h"

// Moose Includes
#include "RayTracing.h"
#include "Ray.h"
#include "MooseError.h"
#include "FEProblem.h"
#include "MooseMesh.h"

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

namespace TraceRay
{

unsigned long int ray_count = 0;

unsigned int debug_ray = 0;

unsigned int debug_ray_id = 0;

unsigned int debug_ray_pid = 0;

// #define USE_DEBUG_RAY

/**
 * Derived from: http://stackoverflow.com/a/565282/2042320
 * q -> q + u*s: Ray
 * p -> p + t*r: Side
 */
template <typename T>
int
sideIntersectedByLine2D(const Elem * current_elem,
                        unsigned int incoming_side,
                        const Point & incoming_point,
                        const std::shared_ptr<Ray> & ray,
                        Point & intersection_point,
                        Point & /*boundary_intersection_point*/)
{
  int intersected_side = -1;
  unsigned int n_sides = current_elem->n_sides();

  std::vector<double> stuff;

  {
//    int boundary_side = -1;

#ifdef USE_DEBUG_RAY
    if (ray->id() == debug_ray_id)
    {
      std::cerr << "Incoming point: " << incoming_point << std::endl;
      std::cerr << "Ray end: " << ray->end() << std::endl;
    }

#endif

    // Bump the starting point down the path a bit
    double q0 = incoming_point(0) + 1e-9 * (ray->end()(0) - incoming_point(0));
    double q1 = incoming_point(1) + 1e-9 * (ray->end()(1) - incoming_point(1));

    double s0 = ray->end()(0) - q0;
    double s1 = ray->end()(1) - q1;

    double intersection_distance = std::numeric_limits<Real>::max();
    double centroid_distance = std::numeric_limits<Real>::max();

    /*
    if (ray_count == 2071730)
      libMesh::err<<std::endl;
    */

    for (unsigned int i = 0; i < n_sides; i++)
    {
      if (i == incoming_side) // Don't search backwards
        continue;

/*
//if (ray_count == 2071730)
{
  libMesh::err<<"side: "<<i<<std::endl;
}
*/

//      bool intersect = false;

#ifdef USE_DEBUG_RAY
      if (ray->id() == debug_ray_id)
      {
        std::cerr << "Side nodes: " << std::endl;
        std::cerr << " " << current_elem->point(T::side_nodes_map[i][0]) << std::endl;
        std::cerr << " " << current_elem->point(T::side_nodes_map[i][1]) << std::endl;
      }
#endif

      double p0 = current_elem->point(T::side_nodes_map[i][0])(0);
      double p1 = current_elem->point(T::side_nodes_map[i][0])(1);

      double r0 = current_elem->point(T::side_nodes_map[i][1])(0) - p0;
      double r1 = current_elem->point(T::side_nodes_map[i][1])(1) - p1;

      double rxs = r0 * s1 - r1 * s0;
      double sxr = s0 * r1 - s1 * r0;

      if (std::abs(rxs) < 1e-10 || std::abs(sxr) < 1e-10) // Lines are parallel or colinear
        continue;

      double t = (((q0 - p0) * s1) - ((q1 - p1) * s0)) / rxs;
      double u = (((p0 - q0) * r1) - ((p1 - q1) * r0)) / sxr;

      /*
      //      if (ray_count == 2071730)
            {
              libMesh::err<<"Tolerance: "<<TOLERANCE<<std::endl;

              libMesh::err<<"t: "<<t<<std::endl;
              libMesh::err<<"u: "<<u<<std::endl;

              libMesh::err<<"t + TOLERANCE: "<<t + TOLERANCE<<std::endl;
              libMesh::err<<"t - TOLERANCE: "<<t - TOLERANCE<<std::endl;
              libMesh::err<<"u + TOLERANCE: "<<u + TOLERANCE<<std::endl;
              libMesh::err<<"u - TOLERANCE: "<<u - TOLERANCE<<std::endl;
            }
      */

      if ((0 < t + 4e-9 && t - 4e-9 <= 1.0) &&
          (0 < u + 4e-9 && u - 4e-9 <= 1.0)) // Do they intersect
      {
        double current_centroid_distance = std::abs(
            t - 0.5); // We want to prefer intersections that go through the middle of sides

        if (current_centroid_distance < centroid_distance)
        {
          //          Elem * neighbor = current_elem->neighbor(i);

          /*
          //if (ray_count == 2071730)
            libMesh::err<<"Intersected side: "<<i<<std::endl;
          */

          intersected_side = i;
          intersection_distance = u;
          centroid_distance = current_centroid_distance;

          // Bump just a bit down the Ray to reduce edge cases
          intersection_point(0) = q0 + ((u /*+1e-9*/) * s0);
          intersection_point(1) = q1 + ((u /*+1e-9*/) * s1);
        }
      }
      else
      {
#ifdef USE_DEBUG_RAY
        if (ray->id() == debug_ray_id)
        {
          std::cerr << "\nSide: " << i << std::endl;
          std::cerr << " 0 < t + 4e-9: " << t + 4e-9 << std::endl;
          std::cerr << " t - 4e-9 <= 1.0: " << t - 4e-9 << std::endl;
          std::cerr << " 0 < u + 4e-9: " << u + 4e-9 << std::endl;
          std::cerr << " u - 4e-9 <= 1.0: " << u - 4e-9 << std::endl;
        }
#endif
      }
    }
  }

  return intersected_side;
}

// https://people.cs.kuleuven.be/~ares.lagae/publications/LD05ERQIT/LD05ERQIT_code.cpp
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
    if (ray->id() == debug_ray_id)
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
    if (ray->id() == debug_ray_id)
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
    if (ray->id() == debug_ray_id)
      libMesh::err << "3 Rejecting because " << t << " < " << 0 << std::endl;
#endif

    return false;
  }

  auto beta = (D * Q) / det;

  if (beta < -1e-12)
  {
#ifdef USE_DEBUG_RAY
    if (ray->id() == debug_ray_id)
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
      if (ray->id() == debug_ray_id)
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
      if (ray->id() == debug_ray_id)
        libMesh::err << "6 Rejecting because " << alpha_prime << " < " << 0 << std::endl;
#endif

      return false;
    }

    auto Q_prime = T_prime.cross(E23);

    auto beta_prime = (D * Q_prime) / det_prime;

    if (beta_prime < -1e-12)
    {
#ifdef USE_DEBUG_RAY
      if (ray->id() == debug_ray_id)
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
                          unsigned int incoming_side,
                          const Point & incoming_point,
                          const std::shared_ptr<Ray> & ray,
                          Point & intersection_point,
                          Point & /*boundary_intersection_point*/)
{

  int intersected_side = -1;
  unsigned int n_sides = 6;

  Point ray_direction = (ray->end() - ray->start());
  ray_direction /= ray_direction.norm();

#ifdef USE_DEBUG_RAY
  if (ray->id() == debug_ray_id)
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
    Point bumped_incoming_point = incoming_point + 1e-9 * ray_direction;

#ifdef USE_DEBUG_RAY
    if (ray->id() == debug_ray_id)
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
    Point centroid(0.5, 0.5, 0);

    for (unsigned int i = 0; i < n_sides; i++)
    {
      if (i == incoming_side) // Don't search backwards
        continue;

#ifdef USE_DEBUG_RAY
      if (ray->id() == debug_ray_id)
        libMesh::err << "Trying side: " << i << std::endl;
#endif

      bool intersected = intersectQuad(bumped_incoming_point,
                                       ray_direction,
                                       *current_elem->get_node(Hex8::side_nodes_map[i][0]),
                                       *current_elem->get_node(Hex8::side_nodes_map[i][1]),
                                       *current_elem->get_node(Hex8::side_nodes_map[i][2]),
                                       *current_elem->get_node(Hex8::side_nodes_map[i][3]),
                                       u_v(0),
                                       u_v(1),
                                       t,
                                       ray);

      if (intersected)
      {
        if (t > intersection_distance)
        {
          //          Elem * neighbor = current_elem->neighbor(i);

          intersected_side = i;
          intersection_distance = t;

          intersection_point = incoming_point + t * ray_direction;
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
  if (ray->id() == debug_ray_id)
  {
    libMesh::err << "Found side: " << intersected_side << std::endl;
    libMesh::err << "Intersection point: " << intersection_point << std::endl;
  }
#endif
  return intersected_side;
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
  unsigned int n_sides = elem->n_sides();

  // Whether or not they intersect
  bool intersect = false;

  unsigned int dim = elem->dim();

  for (unsigned int i=0; i<n_sides; i++)
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
int
sideNeighborIsOn(const Elem * elem, const Elem * neighbor)
{
  unsigned int n_sides = elem->n_sides();

  for (unsigned int i = 0; i < n_sides; i++)
  {
    if (elem->neighbor(i) == neighbor)
      return i;
  }

  return -1;
}

void
find_point_neighbors(const Elem * current_elem,
                     const Point & p,
                     std::set<const Elem *> & neighbor_set,
                     const std::shared_ptr<Ray> & /* ray */)
{
  libmesh_assert(current_elem->contains_point(p));
  libmesh_assert(current_elem->active());

  neighbor_set.clear();
  neighbor_set.insert(current_elem);

  std::set<const Elem *> untested_set, next_untested_set;
  untested_set.insert(current_elem);

#ifdef USE_DEBUG_RAY
  if (ray->id() == debug_ray_id)
    std::cerr << "find_point_neighbors(" << current_elem->id() << ")" << std::endl;
#endif

  while (!untested_set.empty())
  {
    // Loop over all the elements in the patch that haven't already
    // been tested
    std::set<const Elem *>::const_iterator it = untested_set.begin();
    const std::set<const Elem *>::const_iterator end = untested_set.end();

    for (; it != end; ++it)
    {
      const Elem * elem = *it;

      for (unsigned int s = 0; s < elem->n_sides(); s++)
      {
        const Elem * current_neighbor = elem->neighbor_ptr(s);
#ifdef USE_DEBUG_RAY
        if (current_neighbor && ray->id() == debug_ray_id)
          std::cerr << " Testing neighbor " << current_neighbor->id() << std::endl;
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
              if (!neighbor_set.count(current_neighbor))
                next_untested_set.insert(current_neighbor);

#ifdef USE_DEBUG_RAY
              if (ray->id() == debug_ray_id)
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
            std::vector<const Elem *> active_neighbor_children;

            current_neighbor->active_family_tree_by_neighbor(active_neighbor_children, elem);

            std::vector<const Elem *>::const_iterator child_it = active_neighbor_children.begin();
            const std::vector<const Elem *>::const_iterator child_end =
                active_neighbor_children.end();
            for (; child_it != child_end; ++child_it)
            {
              const Elem * current_child = *child_it;
              if (current_child->contains_point(p))
              {
                // Make sure we'll test it
                if (!neighbor_set.count(current_child))
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

/**
 */
void
traceRay(const std::shared_ptr<Ray> & ray,
         RayProblemBase & ray_problem,
         const MooseMesh & mesh,
         unsigned int halo_size,
         Real ray_max_distance,
         Real ray_length,
         THREAD_ID tid)
{
  auto & b_box = ray_problem.boundingBox();

  RaySystem & ray_system = ray_problem.raySystem();

  Point incoming_point = ray->start();

  Point intersection_point;
  Point boundary_intersection_point;
  Point work_point;
  Point work_point2;

  const Elem * current_elem = ray->startingElem();

  unsigned int incoming_side = ray->incomingSide();

  //  unsigned int dim = current_elem->dim();
  unsigned int n_sides = current_elem->n_sides();

  bool keep_tracing = true;

  unsigned long int count = 0;

  unsigned int halo_count = 0;

  ray_count++;

  const std::vector<MooseSharedPointer<RayKernel>> * ray_kernels;

  SubdomainID current_subdomain = -1;

  Mesh * debug_mesh = NULL;
  unsigned int debug_node_count = 0;

  std::vector<BoundaryID> ids(1);

  std::set<const Elem *> point_neighbors;

  // The boundary IDs that have been applied
  std::set<BoundaryID> applied_ids;

  auto pid = mesh.comm().rank();

  //#ifdef USE_DEBUG_RAY
  //  if (pid != debug_ray_pid)
  //    sleep(1000);
  //#endif

  auto mesh_dim = mesh.getMesh().mesh_dimension();

#ifdef USE_DEBUG_RAY
  if (true && ray->id() == debug_ray_id && pid == 0)
  //  if (pid == 1)
  {
    debug_mesh = new Mesh(Parallel::Communicator(), mesh.getMesh().mesh_dimension());
    debug_mesh->skip_partitioning(true);
  }
#endif

//  libMesh::err<<std::endl<<"Starting new ray!"<<std::endl;
#ifdef USE_DEBUG_RAY
  // Add an element for the ray

  if (debug_mesh && ray->id() == debug_ray_id && pid == 0)
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
    int intersected_side = -1;

#ifdef USE_DEBUG_RAY
    if (ray->id() == debug_ray_id)
      std::cerr << "Debug ray UID: " << ray->id() << std::endl;
#endif

#ifdef USE_DEBUG_RAY
    if (ray->id() == debug_ray_id)
      std::cerr << "current_elem: " << current_elem->id() << std::endl;
#endif

    bool ends_in_elem = false;

    // If the ray ends within the mesh... let's see if it ends within this element!
    if (ray->endsWithinMesh())
    {
#ifdef USE_DEBUG_RAY
      if (ray->id() == debug_ray_id)
        std::cerr << "ends within mesh!" << std::endl;
#endif

      if (current_elem->contains_point(ray->end()))
      {
#ifdef USE_DEBUG_RAY
        if (ray->id() == debug_ray_id)
          std::cerr << " ends within current elem!" << std::endl;
#endif
        ends_in_elem = true;
        intersection_point = ray->end();
      }
    }

    if (!ends_in_elem)
    {

      /*

        libMesh::err<<"In element: "<<current_elem->id()<<std::endl;

    if (debug_mesh && ray->id() == debug_ray_id && pid == debug_ray_pid)
  //  if ( == 1)
      {
        for (unsigned int i=0; i<current_elem->n_nodes(); i++)
        {
          debug_mesh->add_point(*current_elem->get_node(i));
          debug_node_count++;
        }

  //      Elem * debug_elem = Elem::build(QUAD4).release();
        Elem * debug_elem = Elem::build(HEX8).release();
        debug_elem->subdomain_id() = 2;

        debug_elem = debug_mesh->add_elem(debug_elem);

        for (unsigned int i=0; i<current_elem->n_nodes(); i++)
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
      for (unsigned int s=0; s<current_elem->n_sides(); s++)
      {
        fe_face->reinit(current_elem, s);
        libMesh::err<<"Side "<<s<<" normal: "<<normals[0]<<std::endl;
      }
    }

      */
      if (current_elem->type() == QUAD4)
        intersected_side = sideIntersectedByLine2D<Quad4>(current_elem,
                                                          incoming_side,
                                                          incoming_point,
                                                          ray,
                                                          intersection_point,
                                                          boundary_intersection_point);
      else if (current_elem->type() == TRI3)
        intersected_side = sideIntersectedByLine2D<Tri3>(current_elem,
                                                         incoming_side,
                                                         incoming_point,
                                                         ray,
                                                         intersection_point,
                                                         boundary_intersection_point);
      if (current_elem->type() == HEX8)
        intersected_side = sideIntersectedByLineHex8(current_elem,
                                                     incoming_side,
                                                     incoming_point,
                                                     ray,
                                                     intersection_point,
                                                     boundary_intersection_point);
    }

    if (intersected_side == -1 && !ends_in_elem) // If we failed to find a side... try harder (as
                                                 // long as the ray has moved a bit
    {
#ifdef USE_DEBUG_RAY
      if (ray->id() == debug_ray_id)
        std::cerr << "Didn't find a side... trying harder" << std::endl;
#endif

      if (ray->distance() > 1e-9) // Only allow this if the ray has moved a bit
      {
        // See if we're on the edge of the domain
        work_point = intersection_point - b_box.min();
        work_point2 = intersection_point - b_box.max();

        unsigned int num_zero = 0;

        for (unsigned int i = 0; i < mesh_dim; i++)
          if (MooseUtils::absoluteFuzzyEqual(std::abs(work_point(i)), 0., 1e-8))
            num_zero++;

        for (unsigned int i = 0; i < mesh_dim; i++)
          if (MooseUtils::absoluteFuzzyEqual(std::abs(work_point2(i)), 0., 1e-8))
            num_zero++;

        if (num_zero) // If we are, then let's just call this point good
        {
#ifdef USE_DEBUG_RAY
//
//        libMesh::err<<ray_count<<" On domain boundary!"<<std::endl;
#endif
          //        auto side = -1;

          auto n_sides = current_elem->n_sides();

          for (auto s = 0u; s < n_sides; s++)
          {
            if (s != incoming_side) // Don't allow a ray to turn around!
            {
              if (current_elem->neighbor(s) ==
                  NULL) // Only check sides that are up against the domain edge
              {
                auto side_elem = current_elem->build_side(s);

                if (side_elem->contains_point(incoming_point))
                {
#ifdef USE_DEBUG_RAY
//
//                libMesh::err<<ray_count<<" Found boundary side: "<<s<<std::endl;
//                libMesh::err<<ray_count<<" Boundary intersection point:
//                "<<incoming_point<<std::endl;
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

      if (intersected_side == -1) // Need to do a more exhaustive search
      {
        // Let's first try grabbing the point_neighbors for this element to see if they are good
        // candidates
        find_point_neighbors(current_elem, incoming_point, point_neighbors, ray);

#ifdef USE_DEBUG_RAY
        if (ray->id() == debug_ray_id)
          std::cerr << "Num point_neighbors: " << point_neighbors.size() << std::endl;
#endif

        const Elem * best_neighbor = NULL;
        unsigned int best_side = 0;
        Real longest_distance = 0;

        // Try the search on each neighbor and see if something good happens
        for (auto & neighbor : point_neighbors)
        {
          if (neighbor != current_elem) // Don't need to look through this element again
          {
#ifdef USE_DEBUG_RAY
            if (ray->id() == debug_ray_id)
              std::cerr << "Trying neighbor " << neighbor->id() << std::endl;
#endif

            // First find the side the point is on in the neighbor
            auto side = -1;

            auto n_sides = neighbor->n_sides();

            for (auto s = 0u; s < n_sides; s++)
            {
              auto side_elem = neighbor->build_side(s);

              if (side_elem->contains_point(incoming_point))
              {
                side = s;
                break;
              }
            }

            // Now try to find an intersection in that neighbor
            if (current_elem->type() == QUAD4)
              intersected_side = sideIntersectedByLine2D<Quad4>(neighbor,
                                                                side,
                                                                incoming_point,
                                                                ray,
                                                                intersection_point,
                                                                boundary_intersection_point);
            else if (current_elem->type() == TRI3)
              intersected_side = sideIntersectedByLine2D<Tri3>(neighbor,
                                                               side,
                                                               incoming_point,
                                                               ray,
                                                               intersection_point,
                                                               boundary_intersection_point);
            if (current_elem->type() == HEX8)
              intersected_side = sideIntersectedByLineHex8(neighbor,
                                                           side,
                                                           incoming_point,
                                                           ray,
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

              //              current_elem = neighbor;
              //              incoming_side = side;

              //              break;
            }
          }
        }

        if (best_neighbor)
        {
#ifdef USE_DEBUG_RAY
          if (debug_ray_id == ray->id())
            std::cerr << "best_neighbor: " << best_neighbor->id() << " best_side: " << best_side
                      << std::endl;
#endif
          current_elem = best_neighbor;
          incoming_side = best_side;
        }

        if (current_elem->processor_id() != pid)
        {
          ray->setStartingElem(current_elem);
          ray->setIncomingSide(incoming_side);
          ray->setStart(incoming_point);
          break;
        }

        count++;

        if (count > 200000)
        {

#ifdef USE_DEBUG_RAY
          if (debug_mesh && ray->id() == debug_ray_id && pid == debug_ray_pid)
          //  if (pid == 1)
          {

            debug_mesh->prepare_for_use();
            debug_mesh->write("debug.e");
          }
#endif
          libMesh::err << pid << " "
                       << "ray_count: " << ray_count << " ray_id: " << ray->id()
                       << " Endless loop while trying harder!" << std::endl;
          libMesh::err << "STUUUUFFFF" << std::endl;
          ray->setShouldContinue(false);

          //      std::abort();

          break;
        }

        continue;
      }
    }

    ray->intersections()++;

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

    if (intersected_side != -1 || ends_in_elem) // -1 means that we didn't find any side
    {
#ifdef USE_DEBUG_RAY
      if (ray->id() == debug_ray_id)
        std::cerr << "Found a side!" << std::endl;
#endif

      //      libMesh::err<<"Intersected side "<<intersected_side<<" at
      //      "<<intersection_point<<std::endl;

      // Use a work point here so temporaries don't get created
      work_point = intersection_point;
      work_point -= incoming_point;

      // The 0.79... is for TY quadrature in 2D
      // TODO: fix this better for 3D
      ray->distance() += work_point.norm();
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
        ray_system.subdomainSetup(current_subdomain, tid);
        ray_kernels = &ray_system.getRayKernels(current_subdomain, tid);

        for (auto & ray_kernel : *ray_kernels)
          ray_kernel->setRay(ray);
      }

      // The "true" causes only sigmaT to be reinited
      ray_system.reinitElem(current_elem, tid, true);

      for (auto & ray_kernel : *ray_kernels)
        ray_kernel->onSegment(current_elem, incoming_point, intersection_point, ends_in_elem);

      // Check if a kernel set to kill the ray
      if (!ray->shouldContinue())
        break;

      // If the Ray is beyond the max distance we need to kill it
      if (ray->distance() >= ray_max_distance)
      {
#ifdef USE_DEBUG_RAY
        if (ray->id() == debug_ray_id)
          std::cerr << "Killing because at max distance!" << std::endl;
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

      if (intersected_side != -1)
      {
        // Get the neighbor on that side
        Elem * neighbor = current_elem->neighbor(intersected_side);

#ifdef USE_DEBUG_RAY
        if (ray->id() == debug_ray_id)
          std::cerr << "Neighbor: " << neighbor << std::endl;
#endif

        if (neighbor)
        {
#ifdef USE_DEBUG_RAY
          if (ray->id() == debug_ray_id)
            std::cerr << "Neighbor_id: " << neighbor->id() << std::endl;
#endif
          // Note: This is finding the side the current_elem is on for the neighbor.  That's the
          // "incoming_side" for the neighbor
          incoming_side = sideNeighborIsOn(neighbor, current_elem);

          // Recurse
          current_elem = neighbor;
          incoming_point = intersection_point;

          // If the neighbor is off-processor... we need to set up the Ray data to be passed to
          // the
          // next processor
          if (neighbor->processor_id() != pid)
          {
            halo_count++;

            if (halo_count == halo_size)
            {
#ifdef USE_DEBUG_RAY
              if (ray->id() == debug_ray_id)
                std::cerr << "Ray going off processor!" << std::endl;
#endif
              ray->setStartingElem(neighbor);
              ray->setIncomingSide(incoming_side);
              ray->setStart(incoming_point);
              break;
            }
          }
        }
        else // Edge of the domain
        {
#ifdef USE_DEBUG_RAY
          if (ray->id() == debug_ray_id)
            std::cerr << "At edge of domain!" << std::endl;
#endif

          /*

                    libMesh::err<<"Reflecting!"<<std::endl;
          */

          //        libMesh::err<<"Original:\n Start: "<<ray->start()<<"\n End:
          //        "<<ray->end()<<std::endl;

          // Need to see if this is a reflective BC
          // Get the "direction" of the Ray
          // Use work_point so no temporaries are created

          //        libMesh::err<<" Direction: "<<work_point<<std::endl;

          mesh.getMesh().get_boundary_info().boundary_ids(current_elem, intersected_side, ids);

          for (const auto & bid : ids)
          {
            if (ray_system.hasRayBCs(bid, tid))
            {
              applied_ids.insert(bid);

              const std::vector<MooseSharedPointer<RayBC>> & ray_bcs =
                  ray_system.getRayBCs(bid, tid);

              for (auto & bc : ray_bcs)
              {
#ifdef USE_DEBUG_RAY
                if (ray->id() == debug_ray_id)
                  std::cerr << ray_count << " Applying BC!" << std::endl;
#endif

                bc->apply(current_elem, intersected_side, intersection_point, ray);
              }
            }
          }

          // Need to see if the ray hit _right_ on the "corner" of the domain
          // If it did, then we're going to apply all of the boundary conditions for each
          // side at that corner

          // First determine if we hit a corner.
          // Subtract off the current point from the min/max of the bounding box
          // If any two entries from that subtraction are (near) zero then we're at a corner
          work_point = intersection_point - b_box.min();
          work_point2 = intersection_point - b_box.max();

          unsigned int num_zero = 0;

          for (unsigned int i = 0; i < mesh_dim; i++)
            if (MooseUtils::absoluteFuzzyEqual(std::abs(work_point(i)), 0., 1e-8))
              num_zero++;

          for (unsigned int i = 0; i < mesh_dim; i++)
            if (MooseUtils::absoluteFuzzyEqual(std::abs(work_point2(i)), 0., 1e-8))
              num_zero++;

          //        if (ray->id() == 4811)
          //          libMesh::err<<"Here!"<<std::endl;

          if (num_zero > 1)
          {
            //          if (ray->id() == 4811)
            //          libMesh::err<<ray_count<<" Ray hit a domain corner!
            //          "<<intersection_point<<std::endl;

            find_point_neighbors(current_elem, intersection_point, point_neighbors, ray);

            // Try the search on each neighbor and see if something good happens
            for (auto & neighbor : point_neighbors)
            {
//            if (neighbor != current_elem) // Don't need to look through this element again
//            {
#ifdef USE_DEBUG_RAY
              if (ray->id() == 4811)
                libMesh::err << "Looking at neighbor: " << neighbor->id() << std::endl;
#endif

              // First find the side the point is on in the neighbor
              auto side = -1;

              auto n_sides = neighbor->n_sides();

              for (auto s = 0u; s < n_sides; s++)
              {
                auto side_elem = neighbor->build_side(s);

                if (side_elem->contains_point(intersection_point))
                {
                  side = s;

                  mesh.getMesh().get_boundary_info().boundary_ids(neighbor, side, ids);

                  for (const auto & bid : ids)
                  {
                    if (applied_ids.find(bid) ==
                        applied_ids.end()) // Don't reapply the same BC twice!
                    {
                      if (ray_system.hasRayBCs(bid, tid))
                      {
                        applied_ids.insert(bid);

                        const std::vector<MooseSharedPointer<RayBC>> & ray_bcs =
                            ray_system.getRayBCs(bid, tid);

                        for (auto & bc : ray_bcs)
                        {
                          //                          libMesh::err<<"Applying BC!"<<std::endl;

                          bc->apply(neighbor, side, intersection_point, ray);
                        }
                      }
                    }
                    //                }
                  }
                }
              }
            }
          }

          applied_ids.clear();

          // See if Ray was killed by the BC
          if (!ray->shouldContinue())
            break;
          else // TODO: This isn't quite right - should really read what was set by the BC
               // later...
          // but this will work for reflective for now
          {
            incoming_point = intersection_point;
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
          //        libMesh::err<<" Direction: "<<work_point<<std::endl;

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

      mooseError("Ray ", ray_count, " killed prematurely!");
      std::cerr << "Ray " << ray_count << " killed prematurely! "
                << " Unique ID: " << ray->id() << std::endl;
      std::cerr << "STUUUUFFFF" << std::endl;

#ifdef USE_DEBUG_RAY
      if (debug_mesh && ray->id() == debug_ray_id && pid == debug_ray_pid)
      //  if (pid == 1)
      {

        debug_mesh->prepare_for_use();
        debug_mesh->write("debug.e");
      }
#endif

      std::abort();

      ray->setShouldContinue(false);
      keep_tracing = false;
    }

    count++;

    if (count > 2000000)
    {

#ifdef USE_DEBUG_RAY
      if (debug_mesh && ray->id() == debug_ray_id && pid == debug_ray_pid)
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

      std::abort();

      break;
    }
  } while (keep_tracing);
}
}
