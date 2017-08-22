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
  params.addParam<bool>("ifenergy",false,"whether to calculate energy");
  params.addParam<bool>("ifstress",false,"whether to calculate stress & jacobian");

  return params;
}

BarcelonaLinearIsoElasticPFDamage::BarcelonaLinearIsoElasticPFDamage(const InputParameters & parameters) :
    ComputeStressBase(parameters),
    _ifenergy(getParam<bool>("ifenergy")),
    _ifstress(getParam<bool>("ifstress")),
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
  Real gc = _gc_prop[_qp];
  Real _m = 1.5 * _Emod[_qp] * gc / (_sigmac[_qp]*_sigmac[_qp]*_l);
  Real lambda = _elasticity_tensor[_qp](0,0,1,1);
  Real mu = _elasticity_tensor[_qp](0,1,0,1);


  RankTwoTensor stress0pos, stress0neg, stress0;
  //Isotropic elasticity is assumed
  Real c = _c[_qp];
  Real _a = (1.0-c)*(1.0-c);
  Real _b = 1.0 + (_m-2.0)*c + (1.0+_p*_m)*c*c;
  Real _degrad = _a / _b;
  Real xfac = _degrad + _kdamage;

  stress0 = _elasticity_tensor[_qp] * _mechanical_strain[_qp];
 

  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

  printf("Successful Decomposition\n");

  Real etr = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    etr += _eigval[i];

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      _eigval[i] = lambda*etr + 2*mu*_eigval[i];

  //Tensors of outerproduct of eigen vectors
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        _etens[i](j,k) = _eigvec(j,i) * _eigvec(k,i);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
      if (_eigval[i] >= 0.0)
        stress0pos += _etens[i] * _eigval[i];
      else
        stress0neg += _etens[i] * _eigval[i];
  }

  //Damage associated with positive component of stress
  _stress[_qp] = stress0pos * xfac + stress0neg;

  //Energy with positive principal strains

if (_ifenergy){

  Real G0_drive = 0.5 * stress0pos.doubleContraction(_mechanical_strain[_qp]);


  Real G0_bar = 0.5 * _sigmac[_qp] * _sigmac[_qp] / _Emod[_qp];

  Real G0_trial =  std::max(G0_bar,G0_drive);

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

if (_ifstress){
  //Used in StressDivergencePFFracTensors Jacobian
  int alpha, beta, i, j, k, l;
  int dim = LIBMESH_DIM;

  //Compute the second derivatives of the eigenvalues w.r.t. strain tensor
  std::vector<RankFourTensor> d2eigvals(dim);
  RankFourTensor dstress_plus, dstress_minus;
  RankFourTensor _Jacob_plus,_Jacob_minus;

  for (alpha = 0; alpha < dim; ++alpha)
      for (beta = 0; beta < dim; ++beta)
      {
          if (_eigval[alpha] == _eigval[beta])
              continue;

          for (i = 0; i < dim; ++i)
              for (j = 0; j < dim; ++j)
                  for (k = 0; k < dim; ++k)
                      for (l = 0; l < dim; ++l)
                      {
                          d2eigvals[alpha](i,j,k,l) += 0.5 * ( _eigvec(beta,i) * _eigvec(alpha,j) + _eigvec(alpha,i) * _eigvec(beta,j) )
                                  * ( _eigvec(beta,k) * _eigvec(alpha,l) + _eigvec(beta,l) * _eigvec(alpha,k) )
                                  / ( _eigval[alpha] - _eigval[beta] );
                      }
      }

      for (alpha = 0; alpha < dim; ++alpha)//loop in eigenvalues
      {
          for (i = 0; i < dim; ++i){
              for (j = 0; j < dim; ++j){
                  for (k = 0; k < dim; ++k){
                      for (l = 0; l < dim; ++l){
                          if (_eigval[alpha]>=0.0) dstress_plus(i,j,k,l) += _eigval[alpha]*d2eigvals[alpha](i,j,k,l) + _etens[alpha](i,j)*_etens[alpha](k,l);
                          if (_eigval[alpha]<0.0) dstress_minus(i,j,k,l) += _eigval[alpha]*d2eigvals[alpha](i,j,k,l) + _etens[alpha](i,j)*_etens[alpha](k,l);
                      }
                  }
              }
          }
      }

      _Jacobian_mult[_qp] = ( xfac * dstress_plus + dstress_minus ) * _elasticity_tensor[_qp];
  }

}
