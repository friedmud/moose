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

#ifndef RAYSYSTEM_H
#define RAYSYSTEM_H

// local includes
#include "RayTracingStudy.h"

// MOOSE includes
#include "SystemBase.h"
#include "ExecuteMooseObjectWarehouse.h"
#include "ConsoleStreamInterface.h"

// libMesh include
#include "libmesh/explicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/petsc_vector.h"

// Forward declarations
class RayProblem;
class TimeIntegrator;
class AuxScalarKernel;
class AuxKernel;
class RayKernel;
class RayBC;
class RayMaterial;

// libMesh forward declarations
namespace libMesh
{
template <typename T>
class NumericVector;
}

/**
 * A system that holds auxiliary variables
 *
 */
class RaySystem : public SystemBase, public ConsoleStreamInterface
{
public:
  /**
   * Data that will be accessed and updated by threads
   *
   * The "alignas(64)" guarantees that this structure will be padded out
   * a 64 byte boundary so that it takes up whole cache lines on modern
   * processors.
   *
   * It also guarantees that a std::vector of these objects will create
   * a contiguous chunk of memory that is 64 byte aligned.
   */
  struct alignas(64) ThreadedData
  {
    /// The current RayMaterial
    MooseSharedPointer<RayMaterial> _current_ray_material;

    /// The index into the solution vectors for the current element
    dof_id_type _current_offset;

    /// The current element being operated on
    const Elem * _current_elem;

    /// Storage for the scalar flux
    PetscVector<Number> * _group_solution;
    PetscScalar * _group_solution_values;

  protected:
    /// Because of portability issues with alignas() > 16 bytes I'm also going to pad with
    /// 64 bytes just to REALLY make sure this stuff is spread out!
    /// Note: this is double (instead of Real) on purpose!
    double padding[8];
  };

  RaySystem(RayProblem & subproblem, const std::string & name, unsigned int num_groups);
  virtual ~RaySystem();

  virtual System & system() { return _sys; };

  virtual NumericVector<Number> & solution() { return *_sys.solution; }
  virtual NumericVector<Number> & solutionOld() { return *_sys.old_local_solution; }
  virtual NumericVector<Number> & solutionOlder() { return *_sys.older_local_solution; }

  virtual void init();

  virtual void initialSetup();
  virtual void timestepSetup();
  virtual void subdomainSetup();
  virtual void residualSetup();
  virtual void jacobianSetup();
  virtual void updateActive(THREAD_ID tid);

  void
  addRayKernel(const std::string & mk_name, const std::string & name, InputParameters parameters);
  void addRayBC(const std::string & mbc_name, const std::string & name, InputParameters parameters);
  void
  addRayMaterial(const std::string & mk_name, const std::string & name, InputParameters parameters);

  /**
   * Called during Ray tracing when the subdomain changes
   * Updates the RayMaterial to the correct one for the subdomain
   */
  void subdomainSetup(SubdomainID current_subdomain, THREAD_ID tid);

  /**
   * Get the RayKernels associated with this subdomain
   */
  const std::vector<MooseSharedPointer<RayKernel>> & getRayKernels(SubdomainID id,
                                                                   THREAD_ID tid = 0) const
  {
    return _ray_kernels.getActiveBlockObjects(id, tid);
  }

  /**
   * Whether or not there are RayBCs for this boundary
   */
  bool hasRayBCs(BoundaryID id, THREAD_ID tid = 0)
  {
    return _ray_bcs.hasActiveBoundaryObjects(id, tid);
  }

  /**
   * Get the RayBCs for this boundary
   */
  const std::vector<MooseSharedPointer<RayBC>> & getRayBCs(BoundaryID id, THREAD_ID tid = 0) const
  {
    return _ray_bcs.getActiveBoundaryObjects(id, tid);
  }

  /**
   * Get the current RayMaterial.
   * Make sure to store the return value as a reference!
   * This is updated automatically when subdomainChanged() is called.
   */
  MooseSharedPointer<RayMaterial> & currentRayMaterial(THREAD_ID tid)
  {
    return _threaded_data[tid]._current_ray_material;
  }

  virtual const NumericVector<Number> *& currentSolution()
  {
    _current_solution = _sys.current_local_solution.get();
    return _current_solution;
  }

  virtual NumericVector<Number> * solutionPreviousNewton() { mooseError("No"); }

  virtual void serializeSolution() { mooseError("No"); }
  virtual NumericVector<Number> & serializedSolution() { mooseError("No"); }

  // This is an empty function since the Aux system doesn't have a matrix!
  virtual void augmentSparsity(SparsityPattern::Graph & /*sparsity*/,
                               std::vector<dof_id_type> & /*n_nz*/,
                               std::vector<dof_id_type> & /*n_oz*/);

  virtual void solve();

  /**
   * Perform one transport sweep
   */
  void transportSweep();

  /**
   * Reinit an element
   * @param elem Which element we are reinitializing for
   * @param tid ID of the thread
   */
  virtual void reinitElem(const Elem * elem, THREAD_ID tid, bool only_sigma_t = false);

  /**
   * The current offset into the solution vectors
   */
  dof_id_type & currentOffset(THREAD_ID tid) { return _threaded_data[tid]._current_offset; }

  /**
   * Get the raw PETSc vector that holds group values for one thread
   * Use this to actually set the value of the scalar flux within RayKernels
   */
  PetscScalar *& groupSolutionValues(THREAD_ID tid)
  {
    return _threaded_data[tid]._group_solution_values;
  }

  /**
   * Get the ray tracing results
   */
  std::map<std::string, std::vector<Real>> & rayTracingResults() { return _results; }

  /**
   * Total number of integrations done in the last transport sweep
   */
  unsigned int totalIntegrations() { return _total_integrations; }

  /**
   * The element currently being operated on by each thread
   */
  const Elem *& currentElem(THREAD_ID tid) { return _threaded_data[tid]._current_elem; }

  /**
   * Set the RayTracingStudy
   */
  void setRayTracingStudy(RayTracingStudy * study) { _study = study; }

protected:
  /**
   * Sweep over all the Rays, applying the Kernels/BCs
   */
  virtual void sweep();

  /**
   * Called after a sweep() has been completed
   *
   * The default is that the values from the threaded group solutions will be accumulated into the
   * true solution vector
   */
  virtual void postSweep();

  RayProblem & _ray_problem;

  TransientExplicitSystem & _sys;

  std::vector<ThreadedData> _threaded_data;

  /// Number of groups
  unsigned int _num_groups;

  /// These are the values that are used internally in the solver
  PetscVector<Number> & _current_group_solution;
  PetscScalar * _current_group_solution_values = NULL;

  /// Storage for RayKernels
  MooseObjectWarehouse<RayKernel> _ray_kernels;
  MooseObjectWarehouse<RayBC> _ray_bcs;
  MooseObjectWarehouse<RayMaterial> _ray_materials;

  /// Ghosted current solution
  const NumericVector<Number> * _current_solution;

  /// Ray tracing information
  std::map<std::string, std::vector<Real>> _results;

  /// The UserObject that propogates Rays
  RayTracingStudy * _study = nullptr;

  /// The amount of integrated distance from the last sweep
  Real _total_integrated_distance;

  /// The total number of "integrations" done during a sweep
  Real _total_integrations;
};

#endif /* EXPLICITSYSTEM_H */
