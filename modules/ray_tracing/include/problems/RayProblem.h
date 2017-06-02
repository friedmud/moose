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
#include "RayAuxSystem.h"
#include "PKESystem.h"

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

  void createRaySystem()
  {
    if (_pke)
      _ray_system = std::make_shared<PKESystem>(*this, "Ray", _num_groups);
    else
      _ray_system = std::make_shared<RaySystem>(*this, "Ray", _num_groups);
  }

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
   * Get the RayAuxSystem for this Problem
   */
  RayAuxSystem & rayAuxSystem() { return _ray_aux_system; }

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
   * Whether or not we're doing TRRM
   */
  bool & isTRRM() { return _trrm; }

  /**
   * Whether or not we're doing a volume only calculation
   */
  bool volumeOnly() { return _volume_only; }

  /**
   * Maximum length through the domain
   */
  Real domainMaxLength() { return _domain_max_length; }

  /**
   * Get the RayTracingStudy for TRRM
   */
  RayTracingStudy & TRRMRayTracingStudy()
  {
    return const_cast<RayTracingStudy &>(
        getUserObject<RayTracingStudy>(_pars.get<UserObjectName>("trrm_study")));
  }

  /**
   * Get the RayTracingStudy
   */
  RayTracingStudy & deterministicRayTracingStudy()
  {
    return const_cast<RayTracingStudy &>(
        getUserObject<RayTracingStudy>(_pars.get<UserObjectName>("deterministic_study")));
  }

  /**
   * Swap the Ray banks
   */
  void swapRayBanks() { std::swap(_current_ray_banks, _old_ray_banks); }

  /**
   * Get the Ray banks
   */
  std::vector<std::vector<std::shared_ptr<Ray>>> & oldRayBanks() { return *_old_ray_banks; }

  /**
   * Get the Ray banks
   */
  std::vector<std::vector<std::shared_ptr<Ray>>> & currentRayBanks() { return *_current_ray_banks; }

  /**
   * Bank a Ray
   */
  void bankRay(const std::shared_ptr<Ray> & ray, THREAD_ID tid)
  {
    (*_current_ray_banks)[tid].emplace_back(ray);
  }

  /**
   * Get a reference to the dead zone value.
   *
   * Note: This is a reference and it should be captured as one since it can change during the
   * simulation.
   */
  const Real & deadZone() { return _dead_zone; }

  /**
   * Get the bounding box for the domain
   */
  const MeshTools::BoundingBox & boundingBox() { return _b_box; }

protected:
  unsigned int _num_groups;
  unsigned int _num_polar;

  bool _trrm;

  bool _deterministic_active;

  bool _volume_only;

  unsigned int _max_rays;

  unsigned int _k_window_size;
  Real _k_window_tolerance;
  unsigned int _k_window_bumps;
  Real _k_tolerance;

  Real _ray_growth;
  Real _ray_length_growth;

  Real _growth_convergence_multiplier;

  Real _scalar_flux_tolerance;

  Real _kappa_source_tolerance;
  Real _kappa_source_growth_tolerance;

  bool _solve_ray;

  bool _solve_fe;
  Real _max_fe_norm;
  Real _fe_tolerance;
  Real _fe_active_tolerance;

  unsigned int _power_iterations;
  unsigned int _conditioning_iterations;

  unsigned int _physics_iterations;

  bool _active_after_converged;
  bool _terminate_after_ray_converged;

  bool _fe_converged = false;
  bool _fe_active = false;
  bool _ray_converged = false;

  Real _dead_zone;

  /// The System that holds all of data for Ray
  std::shared_ptr<RaySystem> _ray_system;

  RayAuxSystem _ray_aux_system;

  Real _domain_max_length;

  /// One for each thread
  std::shared_ptr<std::vector<std::vector<std::shared_ptr<Ray>>>> _current_ray_banks;

  std::shared_ptr<std::vector<std::vector<std::shared_ptr<Ray>>>> _old_ray_banks;

  /// Bounding box for the domain
  MeshTools::BoundingBox _b_box;

  /// Whether or not we're doing a PKE solve
  bool _pke;

  friend class RaySystem;
  friend class RayAuxSystem;
};

#endif /* RAYPROBLEM_H */
