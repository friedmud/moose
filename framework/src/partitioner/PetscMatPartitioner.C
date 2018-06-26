//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PetscMatPartitioner.h"

#include "GeneratedMesh.h"
#include "MooseApp.h"

#include "libmesh/mesh_tools.h"

registerMooseObject("MooseApp", PetscMatPartitioner);

#include <memory>

template <>
InputParameters
validParams<PetscMatPartitioner>()
{
  InputParameters params = validParams<MoosePartitioner>();

  MooseEnum partPackage("parmetis ptscotch chaco party", "ptscotch", false);

  params.addParam<MooseEnum>("part_package",
                             partPackage,
                             "The external package is used for partitioning the mesh via PETSc");

  params.addClassDescription(
      "Partition mesh using external packages via PETSc MatPartitioning interface");

  return params;
}

PetscMatPartitioner::PetscMatPartitioner(const InputParameters & params)
  : MoosePartitioner(params), _part_package(params.get<MooseEnum>("part_package"))
{
}

PetscMatPartitioner::~PetscMatPartitioner() {}

std::unique_ptr<Partitioner>
PetscMatPartitioner::clone() const
{
  return libmesh_make_unique<PetscMatPartitioner>(_pars);
}

void
PetscMatPartitioner::_do_partition(MeshBase & mesh, const unsigned int n_parts)
{
#ifdef LIBMESH_HAVE_PETSC
  std::cout << "Partitioning Mesh!" << std::endl;

  // construct a dual graph
  Mat dual;
  PetscInt *i, *j, *values, *elem_weights, nrows, nj, ncols;
  const PetscInt * parts;
  MatPartitioning part;
  IS is;

  i = j = values = 0;

  // Use libmesh to build the dual graph
  build_graph(mesh);

  nrows = _dual_graph.size();

  // Allocate i: i holds the offset into the adjancency array for each element
  PetscCalloc1(nrows + 1, &i);
  PetscCalloc1(nrows, &elem_weights);

  // Moving index for adjacency array
  nj = 0;


  // Loop over all of the elements and find the _smallest_ side area
  Real min_side_area = std::numeric_limits<Real>::max();
  Real min_elem_surface_area = std::numeric_limits<Real>::max();

  for (auto & elem : mesh.active_local_element_ptr_range())
  {
    Real elem_surface_area = 0;

    for (auto s : elem->side_index_range())
    {
      auto side_area = elem->build_side(s)->volume();

      min_side_area = std::min(side_area, min_side_area);

      elem_surface_area += side_area;
    }

    min_elem_surface_area = std::min(elem_surface_area, min_elem_surface_area);
  }

  // Find the min over all procs
  _communicator.min(min_side_area);
  _communicator.min(min_elem_surface_area);

  // Fill i in with offsets (note that the zeroth entry is zero)
  for (PetscInt k = 0; k < nrows; k++)
  {
    i[k + 1] = i[k] + _dual_graph[k].size();

    auto elem = _local_id_to_elem[k];

    Real elem_surface_area = 0;

    for (auto s : elem->side_index_range())
      elem_surface_area += elem->build_side(s)->volume();

    elem_weights[k] = (elem_surface_area / min_elem_surface_area);
  }

  // j is adjacency matrix
  // The size is the total number of adjancent elements for all elements
  PetscCalloc1(i[nrows], &j);

  // values are edge weights for each adjancency
  PetscCalloc1(i[nrows], &values);

  dof_id_type first_local_elem = 0;
  for (processor_id_type pid = 0; pid < mesh.processor_id(); pid++)
    first_local_elem += _n_active_elem_on_proc[pid];

  {
    unsigned int local_elem_id = 0;

    // Fill in adjancency matrix
    for (auto & row : _dual_graph)
    {
      unsigned int neighbor = 0;

      for (auto adj : row)
      {
        j[nj] = adj;

        auto elem = _local_id_to_elem[local_elem_id];

        auto side_elem = elem->build_side(neighbor);

        values[nj] = (side_elem->volume() / min_side_area);

        nj++;

        neighbor++;
      }

      local_elem_id++;
    }
  }

  ncols = 0;
  for (processor_id_type pid = 0; pid < mesh.n_processors(); pid++)
    ncols += _n_active_elem_on_proc[pid];

  MatCreateMPIAdj(mesh.comm().get(), nrows, ncols, i, j, values, &dual);
  MatPartitioningCreate(mesh.comm().get(), &part);
  MatPartitioningSetAdjacency(part, dual);
  MatPartitioningSetVertexWeights(part, elem_weights);
  MatPartitioningSetNParts(part, n_parts);
  MatPartitioningSetType(part, _part_package.c_str());
  MatPartitioningSetFromOptions(part);
  MatPartitioningApply(part, &is);

  ISGetIndices(is, &parts);

  std::vector<dof_id_type> libmesh_parts;
  std::copy(parts, parts + nrows, std::back_inserter(libmesh_parts));

  ISRestoreIndices(is, &parts);

  assign_partitioning(mesh, libmesh_parts);

  ISRestoreIndices(is, &parts);

  MatPartitioningDestroy(&part);
  ISDestroy(&is);
  MatDestroy(&dual);

  std::cout << "Finished Partitioning Mesh!" << std::endl;

#else
  mooseError("Petsc is required for this partitioner");
#endif
}
