#include "Ray.h"

#include "RayProblem.h"

bool
Ray::operator==(const Ray & other)
{
  return _data == other._data && _start == other._start && _end == other._end &&
         _starting_elem == other._starting_elem && _incoming_side == other._incoming_side &&
         _processor_crossings == other._processor_crossings &&
         _intersections == other._intersections && _distance == other._distance &&
         _integrated_distance == other._integrated_distance &&
         _azimuthal_spacing == other._azimuthal_spacing &&
         _azimuthal_weight == other._azimuthal_weight && _polar_spacing == other._polar_spacing &&
         _polar_sins == other._polar_sins && _polar_weights == other._polar_weights &&
         _ends_within_mesh == other._ends_within_mesh &&
         _ending_elem_id == other._ending_elem_id &&
         _should_continue == other._should_continue;
}

namespace libMesh
{
namespace Parallel
{

unsigned int
Packing<std::shared_ptr<Ray>>::packed_size(typename std::vector<Real>::const_iterator in)
{
  unsigned int total_size = 0;

  // Spot for the data length
  total_size += 1;

  // Spot for the polar_sin length
  total_size += 1;

  // Spot for the polar_weight length
  total_size += 1;

  // Spots for the data
  total_size += *in++;

  // Spots for the polar_sins
  total_size += *in++;

  // Spots for the polar_weights
  total_size += *in++;

  // Points
  total_size += LIBMESH_DIM * 2;

  // Starting element, ending element id, ends within mesh, incoming side,
  // processor crosssings, intersections, distance, integrated_distance,
  // azimuthal_spacing, azimuthal_weight, polar_spacing, id, is_reverse,
  // should continue
  total_size += 14;

  return total_size;
}

unsigned int
Packing<std::shared_ptr<Ray>>::packable_size(const std::shared_ptr<Ray> & ray, const void *)
{
  unsigned int total_size = 0;

  // Spot for the data length
  total_size += 1;

  // Spot for the polar_sin length
  total_size += 1;

  // Spot for the polar_weight length
  total_size += 1;

  // Spots for the data
  total_size += ray->_data.size();

  // Spots for the polar sins
  total_size += ray->_polar_sins.size();

  // Spots for the polar weights
  total_size += ray->_polar_weights.size();

  // Points
  total_size += LIBMESH_DIM * 2;

  // Starting element, ending element id, ends within mesh, incoming side,
  // processor crosssings, intersections, distance, integrated_distance,
  // azimuthal_spacing, azimuthal_weight, polar_spacing, id, is_reverse,
  // should continue
  total_size += 14;

  return total_size;
}

std::shared_ptr<Ray>
Packing<std::shared_ptr<Ray>>::unpack(typename std::vector<Real>::const_iterator in,
                                      void * ray_problem_ptr)
{
  //    auto unpacking_start_time = std::chrono::steady_clock::now();

  auto ray_problem = static_cast<RayProblemBase *>(ray_problem_ptr);

  std::shared_ptr<Ray> ray = ray_problem->_ray_pool.acquire();

  // Grab the data size
  unsigned int data_size = *in++;

  // Grab the polar_sins size
  unsigned int polar_sins_size = *in++;

  // Grab the polar_weights size
  unsigned int polar_weights_size = *in++;

  // Start Point
  ray->_start(0) = (*in++);
  ray->_start(1) = (*in++);
  ray->_start(2) = (*in++);

  // End Point
  ray->_end(0) = (*in++);
  ray->_end(1) = (*in++);
  ray->_end(2) = (*in++);

  // Starting Element
  auto mesh = static_cast<MeshBase *>(&ray_problem->mesh().getMesh());
  ray->setStartingElem(mesh->elem((*in++)));

  // Ending element
  auto id = (*in++);
  if (id == (Real)DofObject::invalid_id)
    ray->setEndingElemId(DofObject::invalid_id);
  else
    ray->setEndingElemId(id);

  // Ends Within Mesh
  ray->setEndsWithinMesh((*in++));

  // Incoming side
  ray->setIncomingSide((*in++));

  // Processor Crossings
  ray->processorCrossings() = (*in++);

  // Intersections
  ray->intersections() = (*in++);

  // Distance
  ray->distance() = (*in++);

  // Integrated Distance
  ray->integratedDistance() = (*in++);

  // Azimuthal Spacing
  ray->setAzimuthalSpacing((*in++));

  // Azimuthal Weight
  ray->setAzimuthalWeight((*in++));

  // Polar Spacing
  ray->setPolarSpacing((*in++));

  // ID
  ray->setID((*in++));

  // is_reverse
  ray->setIsReverse((*in++));

  // Should continue
  ray->setShouldContinue((*in++));

  // Reserve space for the data
  ray->_data.resize(data_size);

  // Copy out data
  for (unsigned int i = 0; i < data_size; i++)
    ray->_data[i] = *in++;

  // Reserve space for the polar_sins
  ray->_polar_sins.resize(polar_sins_size);

  // Copy out the polar_sins
  for (unsigned int i = 0; i < polar_sins_size; i++)
    ray->_polar_sins[i] = *in++;

  // Reserve space for the polar_weights
  ray->_polar_weights.resize(polar_weights_size);

  // Copy out the polar_weights
  for (unsigned int i = 0; i < polar_weights_size; i++)
    ray->_polar_weights[i] = *in++;

  return ray;
}

} // namespace Parallel

} // namespace libMesh
