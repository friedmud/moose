//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterActions.h"

#include "Factory.h"

#include "Syntax.h"
#include "ActionFactory.h"

#include "AddBCAction.h"
#include "AddDiracKernelAction.h"
#include "AddICAction.h"
#include "AddInitialConditionAction.h"
#include "AddKernelAction.h"
#include "AddScalarKernelAction.h"
#include "AddDGKernelAction.h"
#include "AddInterfaceKernelAction.h"
#include "AddPeriodicBCAction.h"
#include "AddVariableAction.h"
#include "AddAuxVariableAction.h"
#include "AddPostprocessorAction.h"
#include "AddVectorPostprocessorAction.h"
#include "AddDamperAction.h"
#include "AddFunctionAction.h"
#include "AddDistributionAction.h"
#include "AddSamplerAction.h"
#include "CreateExecutionerAction.h"
#include "DetermineSystemType.h"
#include "EmptyAction.h"
#include "InitProblemAction.h"
#include "CopyNodalVarsAction.h"
#include "SetupMeshAction.h"
#include "AddMeshModifierAction.h"
#include "SetupMeshCompleteAction.h"
#include "AddOutputAction.h"
#include "CommonOutputAction.h"
#include "AddMaterialAction.h"
#include "GlobalParamsAction.h"
#include "AdaptivityAction.h"
#include "PartitionerAction.h"
#include "SetupDampersAction.h"
#include "CheckIntegrityAction.h"
#include "SetupQuadratureAction.h"
#include "SetupPreconditionerAction.h"
#include "SetupDebugAction.h"
#include "SetupResidualDebugAction.h"
#include "DeprecatedBlockAction.h"
#include "AddConstraintAction.h"
#include "CreateDisplacedProblemAction.h"
#include "CreateProblemAction.h"
#include "DynamicObjectRegistrationAction.h"
#include "AddUserObjectAction.h"
#include "AddControlAction.h"
#include "AddElementalFieldAction.h"
#include "AddIndicatorAction.h"
#include "AddMarkerAction.h"
#include "SetAdaptivityOptionsAction.h"
#include "AddMultiAppAction.h"
#include "AddTransferAction.h"
#include "AddNodalNormalsAction.h"
#include "SetupTimeStepperAction.h"
#include "SetupTimeIntegratorAction.h"
#include "SetupPredictorAction.h"
#include "AddMortarInterfaceAction.h"
#include "SetupPostprocessorDataAction.h"
#include "MaterialOutputAction.h"
#include "CheckOutputAction.h"
#include "SetupRecoverFileBaseAction.h"
#include "AddNodalKernelAction.h"
#include "MaterialDerivativeTestAction.h"
#include "AddRelationshipManager.h"
#include "MeshOnlyAction.h"
#include "SplitMeshAction.h"
#include "AddFieldSplitAction.h"
#include "AddBoundsVectorsAction.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerActionObjects(Syntax & syntax, ActionFactory & action_factory)
{
#undef registerAction
#define registerAction(tplt, action)                                                               \
  action_factory.reg<tplt>(stringifyName(tplt), action, __FILE__, __LINE__)

  registerAction(SetupPostprocessorDataAction, "setup_postprocessor_data");

  registerAction(SplitMeshAction, "split_mesh");
  registerAction(MeshOnlyAction, "mesh_only");
  registerAction(SetupMeshAction, "setup_mesh");
  registerAction(SetupMeshAction, "init_mesh");
  registerAction(SetupMeshCompleteAction, "prepare_mesh");
  registerAction(AddMeshModifierAction, "add_mesh_modifier");
  registerAction(AddMortarInterfaceAction, "add_mortar_interface");
  registerAction(SetupMeshCompleteAction, "execute_mesh_modifiers");
  registerAction(SetupMeshCompleteAction, "uniform_refine_mesh");
  registerAction(SetupMeshCompleteAction, "setup_mesh_complete");

  registerAction(AddFunctionAction, "add_function");
  registerAction(AddDistributionAction, "add_distribution");
  registerAction(AddSamplerAction, "add_sampler");
  registerAction(CreateExecutionerAction, "setup_executioner");
  registerAction(SetupTimeStepperAction, "setup_time_stepper");
  registerAction(SetupTimeIntegratorAction, "setup_time_integrator");
  registerAction(CreateDisplacedProblemAction, "init_displaced_problem");
  registerAction(DetermineSystemType, "determine_system_type");
  registerAction(CreateProblemAction, "create_problem");
  registerAction(DynamicObjectRegistrationAction, "dynamic_object_registration");
  registerAction(AddOutputAction, "add_output");
  registerAction(CommonOutputAction, "common_output");
  registerAction(SetupRecoverFileBaseAction, "setup_recover_file_base");
  registerAction(GlobalParamsAction, "set_global_params");
  registerAction(SetupPredictorAction, "setup_predictor");
  registerAction(MaterialOutputAction, "setup_material_output");
  registerAction(CheckOutputAction, "check_output");

  /// Variable/AuxVariable Actions
  registerAction(AddVariableAction, "add_variable");
  registerAction(AddAuxVariableAction, "add_aux_variable");

  registerAction(CopyNodalVarsAction, "check_copy_nodal_vars");
  registerAction(CopyNodalVarsAction, "copy_nodal_vars");
  registerAction(CopyNodalVarsAction, "copy_nodal_aux_vars");

  // Initial Condition Actions
  registerAction(AddICAction, "add_ic");
  registerAction(AddInitialConditionAction, "add_ic");

  registerAction(AddKernelAction, "add_kernel");
  registerAction(AddNodalKernelAction, "add_nodal_kernel");
  registerAction(AddKernelAction, "add_aux_kernel");
  registerAction(AddScalarKernelAction, "add_scalar_kernel");
  registerAction(AddScalarKernelAction, "add_aux_scalar_kernel");
  registerAction(AddDGKernelAction, "add_dg_kernel");
  registerAction(AddInterfaceKernelAction, "add_interface_kernel");
  registerAction(AddBCAction, "add_bc");
  registerAction(EmptyAction, "no_action"); // placeholder
  registerAction(AddPeriodicBCAction, "add_periodic_bc");
  registerAction(AddMaterialAction, "add_material");
  registerAction(AddPostprocessorAction, "add_postprocessor");
  registerAction(AddVectorPostprocessorAction, "add_vector_postprocessor");
  registerAction(AddDamperAction, "add_damper");
  registerAction(AddFieldSplitAction, "add_field_split");
  registerAction(SetupPreconditionerAction, "add_preconditioning");
  registerAction(SetupQuadratureAction, "setup_quadrature");
  registerAction(DeprecatedBlockAction, "deprecated_block");
  registerAction(AddConstraintAction, "add_constraint");
  registerAction(AddUserObjectAction, "add_user_object");
  registerAction(AddControlAction, "add_control");
  registerAction(AddElementalFieldAction, "add_elemental_field_variable");
  registerAction(AddIndicatorAction, "add_indicator");
  registerAction(AddMarkerAction, "add_marker");
  registerAction(SetAdaptivityOptionsAction, "set_adaptivity_options");

  registerAction(AddNodalNormalsAction, "add_aux_variable");
  registerAction(AddNodalNormalsAction, "add_postprocessor");
  registerAction(AddNodalNormalsAction, "add_user_object");

#ifdef LIBMESH_ENABLE_AMR
  registerAction(AdaptivityAction, "setup_adaptivity");
#endif

  registerAction(PartitionerAction, "add_partitioner");
  registerAction(AddDiracKernelAction, "add_dirac_kernel");
  registerAction(SetupDebugAction, "setup_debug");
  registerAction(SetupResidualDebugAction, "setup_residual_debug");

  registerAction(AddBoundsVectorsAction, "add_bounds_vectors");

  // NonParsedActions
  registerAction(SetupDampersAction, "setup_dampers");
  registerAction(EmptyAction, "ready_to_init");
  registerAction(AddRelationshipManager, "add_algebraic_rm");
  registerAction(AddRelationshipManager, "add_geometric_rm");

  registerAction(InitProblemAction, "init_problem");
  registerAction(CheckIntegrityAction, "check_integrity");

  registerAction(AddMultiAppAction, "add_multi_app");
  registerAction(AddTransferAction, "add_transfer");

  // TODO: Why is this here?
  registerTask("finish_input_file_output", false);
  registerAction(EmptyAction, "finish_input_file_output");

  registerAction(MaterialDerivativeTestAction, "add_variable");
  registerAction(MaterialDerivativeTestAction, "add_kernel");
  registerAction(MaterialDerivativeTestAction, "add_preconditioning");

#undef registerAction
#define registerAction(tplt, action) action_factory.regLegacy<tplt>(stringifyName(tplt), action)
}
}
