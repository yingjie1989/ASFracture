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


//include header of my application
#include "ASFracture.h"

#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "ActionFactory.h"
#include "MooseSyntax.h"


//kernels used from modules
#include "TensorMechanicsApp.h"
#include "PhaseFieldApp.h"
#include "HeatConductionApp.h"

//initial condition
//#include "CrossIC.h"

//auxkernel
#include "NewmarkDispAux.h"

//boundary condition
#include "CoupledNeumannBC.h"
#include "CoupledNeumannVectorBC.h"


//material
#include "DerivativeParsedMaterial.h"
#include "LinearIsoElasticPFDamageModify.h"
#include "CohesiveLinearIsoElasticPFDamage.h"
#include "PFFracRandomBulkRateMaterial.h"
#include "WeibullMaterial.h"


//custom kernel
#include "Diffusion_D.h"
#include "CohesivePFFracBulkRate.h"
#include "PFFracBulkRateModify.h"
#include "PFFracBulkRateAxisymmetric.h"
#include "SourceMonopole.h"
#include "StressDivergenceRZPFFracTensors.h"
#include "MassLumpedReaction.h"
#include "StressDivergenceExpPFFracTensors.h"
#include "StressDivergenceExpTensors.h"
#include "StressDivergenceExplicitTensors.h"
#include "InertialForceExp.h"
#include "PFFracIntVar.h"
#include "PFFracCoupledInterfaceExp.h"
#include "TimeDerivativeExp.h"

//dirac kernel
#include "MonopoleDirac.h"

//aux kernel
#include "ExpAccelAux.h"
#include "ExpVelAux.h"


template<>
InputParameters validParams<ASFracture>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

ASFracture::ASFracture(InputParameters parameters) :
    MooseApp(parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ASFracture::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ASFracture::associateSyntax(_syntax, _action_factory);
}

ASFracture::~ASFracture()
{
}

void
ASFracture::registerObjects(Factory & factory)
{
//register objects from module phase_field
  TensorMechanicsApp::registerObjects(factory);
  PhaseFieldApp::registerObjects(factory);
  HeatConductionApp::registerObjects(factory);


//kernels
registerKernel(Diffusion_D);
registerKernel(PFFracBulkRateModify);
registerKernel(CohesivePFFracBulkRate);
registerKernel(PFFracBulkRateAxisymmetric);
registerKernel(SourceMonopole);
registerKernel(StressDivergenceRZPFFracTensors);
registerKernel(MassLumpedReaction);
registerKernel(StressDivergenceExpPFFracTensors);
registerKernel(StressDivergenceExpTensors);
registerKernel(StressDivergenceExplicitTensors);
registerKernel(InertialForceExp);
//registerKernel(PFFracIntVar);
registerKernel(TimeDerivativeExp);
registerKernel(PFFracCoupledInterfaceExp);



//IC
//registerInitialCondition(FingerIC);

//Dirackernels
registerDiracKernel(MonopoleDirac);


//BC
registerBoundaryCondition(CoupledNeumannBC);
registerBoundaryCondition(CoupledNeumannVectorBC);

//Materials

registerMaterial(LinearIsoElasticPFDamageModify);
registerMaterial(CohesiveLinearIsoElasticPFDamage);
registerMaterial(PFFracRandomBulkRateMaterial);
registerMaterial(WeibullMaterial);


//Auxkernels
registerAux(NewmarkDispAux);
registerAux(ExpAccelAux);
registerAux(ExpVelAux);


}

void
ASFracture::registerApps()
{
  registerApp(ASFracture);
}

void
ASFracture::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
   TensorMechanicsApp::associateSyntax(syntax, action_factory);
   PhaseFieldApp::associateSyntax(syntax, action_factory);
   HeatConductionApp::associateSyntax(syntax, action_factory);

}
