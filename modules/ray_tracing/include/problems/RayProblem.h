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

#ifndef RAYPROBLEM_H
#define RAYPROBLEM_H

// Local Includes
#include "RaySystem.h"

// MOOSE Includes
#include "FEProblem.h"
#include "MooseObjectWarehouse.h"

// libMesh Includes
#include "libmesh/petsc_vector.h"

class RayProblem;

class RayTracingStudy;

template <>
InputParameters validParams<RayProblem>();

/**
 * FEProblem derived class for customization of callbacks. In this instance we only print out
 * something in the c-tor and d-tor, so we know the class was build and used properly.
 */
class RayProblem : public FEProblem
{
public:
  RayProblem(const InputParameters & params);
  virtual ~RayProblem() override;

  void createRaySystem() { _ray_system = std::make_shared<RaySystem>(*this, "Ray", _num_groups); }

  virtual void initialSetup() override;
  virtual void timestepSetup() override;

  virtual void solve() override;

  virtual bool converged() override;

  virtual void subdomainSetup(SubdomainID subdomain, THREAD_ID tid) override;

  /**
   * Run subdomain setup on just the Ray systems
   */
  virtual void RaySubdomainSetup(SubdomainID subdomain, THREAD_ID tid);

  virtual void reinitElem(const Elem * elem, THREAD_ID tid) override;
  virtual void
  reinitElemPhys(const Elem * elem, std::vector<Point> phys_points_in_elem, THREAD_ID tid) override;

  /**
   * Run reinitElem() on just the Ray systems
   */
  virtual void RayReinitElem(const Elem * elem, THREAD_ID tid);

  virtual void prepareMaterials(SubdomainID blk_id, THREAD_ID tid) override;
  virtual void
  reinitMaterials(SubdomainID blk_id, THREAD_ID tid, bool swap_stateful = true) override;

  virtual std::vector<VariableName> getVariableNames() override;

  virtual bool hasVariable(const std::string & var_name) override;

  virtual MooseVariable & getVariable(THREAD_ID tid, const std::string & var_name) override;

  /**
   * Get the RaySystem for this Problem
   */
  RaySystem & raySystem() { return *_ray_system; }

  /**
   * Called during the ray tracing process when the subdomain changes
   */
  void subdomainChanged(SubdomainID current_subdomain, THREAD_ID tid);

  /**
   * The number of energy groups currently being utilized.
   */
  unsigned int numGroups() { return _num_groups; }

  /**
   * The number of polar angles currently being used
   */
  unsigned int numPolar() { return _num_polar; }

  /**
   * Maximum length through the domain
   */
  Real domainMaxLength() { return _domain_max_length; }

  /**
   * Get the RayTracingStudy
   */
  virtual RayTracingStudy & rayTracingStudy()
  {
    return const_cast<RayTracingStudy &>(
        getUserObject<RayTracingStudy>(_pars.get<UserObjectName>("trrm_study")));
  }

  /**
   * Get the bounding box for the domain
   */
  const MeshTools::BoundingBox & boundingBox() { return _b_box; }

protected:
  unsigned int _num_groups;
  unsigned int _num_polar;

  bool _solve_ray;

  bool _solve_fe;

  /// The System that holds all of data for Ray
  std::shared_ptr<RaySystem> _ray_system;

  Real _domain_max_length;

  /// Bounding box for the domain
  MeshTools::BoundingBox _b_box;

  friend class RaySystem;
  friend class RayAuxSystem;
};

#endif /* RAYPROBLEM_H */
