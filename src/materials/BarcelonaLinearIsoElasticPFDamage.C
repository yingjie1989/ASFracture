/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "BarcelonaLinearIsoElasticPFDamage.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<BarcelonaLinearIsoElasticPFDamage>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Phase-field fracture model energy contribution to damage growth-isotropic elasticity and undamaged stress under compressive strain");
  params.addRequiredCoupledVar("c","Order parameter for damage");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addParam<bool>("historyEng",false,"indicator whether to use history strain energy");
  params.addRequiredParam<Real>("l","Interface width");
  params.addRequiredParam<Real>("p","p parameter which influences the cohesive traction separation law");
  params.addRequiredParam<MaterialPropertyName>("gc_prop_var", "Material property name with gc value");
  params.addRequiredParam<MaterialPropertyName>("Emod", "Material property name with Young's Modulus");
  params.addRequiredParam<MaterialPropertyName>("sigmac", "Material property name with strength");

  return params;
}

BarcelonaLinearIsoElasticPFDamage::BarcelonaLinearIsoElasticPFDamage(const InputParameters & parameters) :
    ComputeStressBase(parameters),
    _c(coupledValue("c")),
    _kdamage(getParam<Real>("kdamage")),
    _historyEng(getParam<bool>("historyEng")),
    _gc_prop(getMaterialProperty<Real>("gc_prop_var")),
    _Emod(getMaterialProperty<Real>("Emod")),
    _sigmac(getMaterialProperty<Real>("sigmac")),
    _l(getParam<Real>("l")),
    _p(getParam<Real>("p")),
    _G0_pos(declareProperty<Real>("G0_pos")),
    _G0_pos_old(declarePropertyOld<Real>("G0_pos")),
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _dG0_pos_dstrain(declareProperty<RankTwoTensor>("dG0_pos_dstrain")),
    _etens(LIBMESH_DIM),
    _epos(LIBMESH_DIM),
    _eigval(LIBMESH_DIM)
{

}

void
BarcelonaLinearIsoElasticPFDamage::initQpStatefulProperties()
{
   ComputeStressBase::initQpStatefulProperties();

   //printf("time is %lf\n",_t);

   if (_t > 0)
  	 _G0_pos_old[_qp] = _G0_pos[_qp];
   else{
	    _G0_pos[_qp] = 0.0;
	    _G0_pos_old[_qp] = 0.0;
   }
   _dG0_pos_dstrain[_qp] = _stress[_qp];

   _dstress_dc[_qp] = _stress[_qp] * 0.0;

}


void BarcelonaLinearIsoElasticPFDamage::computeQpStress()
{
  updateVar();
  updateJacobian();
}

void
BarcelonaLinearIsoElasticPFDamage::updateVar()
{
  Real gc = _gc_prop[_qp];
  Real _m = 1.5 * _Emod[_qp] * gc / (_sigmac[_qp]*_sigmac[_qp]*_l);


  RankTwoTensor stress0pos, stress0neg, stress0;
  //Isotropic elasticity is assumed
  Real lambda = _elasticity_tensor[_qp](0,0,1,1);
  Real mu = _elasticity_tensor[_qp](0,1,0,1);
  Real c = _c[_qp];
  Real _a = (1.0-c)*(1.0-c);
  Real _b = 1.0 + (_m-2.0)*c + (1.0+_p*_m)*c*c;
  Real _degrad = _a / _b;
  Real xfac = _degrad*(1.0-_kdamage) + _kdamage;

  stress0 = _elasticity_tensor * _mechanical_strain;

  stress0.symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

  //Tensors of outerproduct of eigen vectors
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        _etens[i](j,k) = _eigvec(j,i) * _eigvec(k,i);


  Real etr = 0.0;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
      stress0pos += _etens[i] * (std::abs(_eigval[i]) + _eigval[i]) / 2.0;
      stress0neg += _etens[i] * (std::abs(_eigval[i]) - _eigval[i]) / 2.0;
      if (_eigval[i] > 0)
        etr += _eigval[i] * _eigval[i];
  }

  etr -= _sigmac[_qp] * _sigmac[_qp];

  //Damage associated with positive component of stress
  _stress[_qp] = stress0pos * xfac - stress0neg;

  //Energy with positive principal strains

  Real G0_bar = 0.5 * _sigmac[_qp] * _sigmac[_qp] / _Emod[_qp];

  Real G0_trial = 0.5 / _Emod[_qp] * std::max(etr,0.0) + G0_bar;

 if (!_historyEng){
      _G0_pos[_qp] = G0_trial;
      _dG0_pos_dstrain[_qp] = stress0pos;
 }else{
      if (G0_trial > _G0_pos_old[_qp]){
	        _G0_pos[_qp] = G0_trial;
  	      _dG0_pos_dstrain[_qp] = stress0pos;
      }else{
	       _G0_pos[_qp] = _G0_pos_old[_qp];
	       _dG0_pos_dstrain[_qp] = stress0pos * 0.0;
      }
 }
  //Used in StressDivergencePFFracTensors Jacobian
  Real _da_dphi = - 2.0 * ( 1.0 - c );
  Real _db_dphi = (_m-2.0) + (1.0+_p*_m)*2.0*c;

  Real _dg_dphi = _da_dphi/_b - _a/(_b*_b) * _db_dphi;
  _dstress_dc[_qp] = stress0pos * _dg_dphi;

}

void
BarcelonaLinearIsoElasticPFDamage::updateJacobian()
{
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
