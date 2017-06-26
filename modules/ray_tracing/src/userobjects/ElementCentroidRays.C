#include "ElementCentroidRays.h"

// MOOSE Includes
#include "MooseMesh.h"

template <>
InputParameters
validParams<ElementCentroidRays>()
{
  InputParameters params = validParams<RayTracingStudy>();

  params.addRequiredParam<std::vector<dof_id_type>>("start_elems",
                                                    "The elements to start the Rays from");
  params.addRequiredParam<std::vector<dof_id_type>>("end_elems",
                                                    "The elements to start the Rays from");
  return params;
}

ElementCentroidRays::ElementCentroidRays(const InputParameters & parameters)
  : RayTracingStudy(parameters),
    _start_elems(getParam<std::vector<dof_id_type>>("start_elems")),
    _end_elems(getParam<std::vector<dof_id_type>>("end_elems"))
{
  if (_start_elems.size() != _end_elems.size())
    mooseError("start_elems length must match end_elems length!");
}

void
ElementCentroidRays::generateRays()
{
  std::vector<std::shared_ptr<Ray>> rays;

  for (auto i = beginIndex(_start_elems); i < _start_elems.size(); i++)
  {
    auto start_elem_id = _start_elems[i];
    auto end_elem_id = _end_elems[i];

    auto start_elem = _mesh.elemPtr(start_elem_id);
    auto end_elem = _mesh.elemPtr(end_elem_id);

    auto ray = std::make_shared<Ray>(
        start_elem->centroid(), end_elem->centroid(), _num_groups, start_elem);

    ray->setEndsWithinMesh();

    rays.emplace_back(ray);
  }

  chunkyTraceAndBuffer(rays);
}
