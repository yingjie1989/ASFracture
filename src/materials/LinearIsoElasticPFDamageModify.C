/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "LinearIsoElasticPFDamageModify.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<LinearIsoElasticPFDamageModify>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Phase-field fracture model energy contribution to damage growth-isotropic elasticity and undamaged stress under compressive strain");
  params.addRequiredCoupledVar("c","Order parameter for damage");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");

  return params;
}

LinearIsoElasticPFDamageModify::LinearIsoElasticPFDamageModify(const InputParameters & parameters) :
    ComputeStressBase(parameters),
    _c(coupledValue("c")),
    _kdamage(getParam<Real>("kdamage")),
    _G0_pos(declareProperty<Real>("G0_pos")),
    _G0_pos_old(declarePropertyOld<Real>("G0_pos")),
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _dG0_pos_dstrain(declareProperty<RankTwoTensor>("dG0_pos_dstrain")),
    _etens(LIBMESH_DIM),
    _epos(LIBMESH_DIM),
    _eigval(LIBMESH_DIM)
{

}

void LinearIsoElasticPFDamageModify::initQpStatefulProperties() 
{
   ComputeStressBase::initQpStatefulProperties();

   //printf("time is %lf\n",_t);

   if(_t > 0)
  	 _G0_pos_old[_qp] = _G0_pos[_qp]; 
   else{
	 _G0_pos[_qp] = 0.0;
	 _G0_pos_old[_qp] = 0.0;
   }
   _dG0_pos_dstrain[_qp] = _stress[_qp];
   _dstress_dc[_qp] = -_stress[_qp] * (2.0 * (1.0 - _c[_qp]));  

}


void LinearIsoElasticPFDamageModify::computeQpStress()
{
  updateVar();
  updateJacobian();
}

void
LinearIsoElasticPFDamageModify::updateVar()
{
  RankTwoTensor stress0pos, stress0neg, stress0;
  //Isotropic elasticity is assumed
  Real lambda = _elasticity_tensor[_qp](0,0,1,1);
  Real mu = _elasticity_tensor[_qp](0,1,0,1);
  Real c = _c[_qp];
  Real xfac = ( Utility::pow<2>(1.0-c) )*(1-_kdamage) + _kdamage;

  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

  //Tensors of outerproduct of eigen vectors
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        _etens[i](j,k) = _eigvec(j,i) * _eigvec(k,i);

  Real etr = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    etr += _eigval[i];

  Real etrpos = (std::abs(etr) + etr) / 2.0;
  Real etrneg = (std::abs(etr) - etr) / 2.0;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    stress0pos += _etens[i] * (lambda * etrpos + 2.0 * mu * (std::abs(_eigval[i]) + _eigval[i]) / 2.0);
    stress0neg += _etens[i] * (lambda * etrneg + 2.0 * mu * (std::abs(_eigval[i]) - _eigval[i]) / 2.0);
  }

  //Damage associated with positive component of stress
  _stress[_qp] = stress0pos * xfac - stress0neg;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    _epos[i] = (std::abs(_eigval[i]) + _eigval[i]) / 2.0;

  Real val = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    val += Utility::pow<2>(_epos[i]);
  val *= mu;

  //Before Update

  //Energy with positive principal strains

  Real G0_trial = lambda * Utility::pow<2>(etrpos) / 2.0 + val;

  //printf("material properties is %lf, %lf\n",_G0_pos[_qp],_G0_pos_old[_qp]);  

  if(G0_trial > _G0_pos_old[_qp]){
	_G0_pos[_qp] = G0_trial;
  	_dG0_pos_dstrain[_qp] = stress0pos;
  }else{
	//printf("irreversibility is enforced\n");
	_G0_pos[_qp] = _G0_pos_old[_qp];
	_dG0_pos_dstrain[_qp] = stress0pos * 0.0;
  }

  //Used in StressDivergencePFFracTensors Jacobian
  _dstress_dc[_qp] = -stress0pos * (2.0 * (1.0 - c));
}

void
LinearIsoElasticPFDamageModify::updateJacobian()
{
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
