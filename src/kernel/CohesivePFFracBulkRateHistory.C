/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CohesivePFFracBulkRateHistory.h"

template<>
InputParameters validParams<CohesivePFFracBulkRateHistory>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Kernel to compute bulk energy contribution to damage order parameter residual equation");
  params.addRequiredParam<Real>("l","Interface width");
   params.addRequiredParam<Real>("p","p parameter which influences the cohesive traction separation law");
  params.addRequiredParam<Real>("visco","Viscosity parameter");
  params.addParam<bool>("ifOld",false,"indicator whether to use old value");
  params.addRequiredParam<MaterialPropertyName>("gc_prop_var", "Material property name with gc value");
  params.addRequiredParam<MaterialPropertyName>("G0_var", "Material property name with undamaged strain energy driving damage (G0_pos)");
  params.addParam<MaterialPropertyName>("dG0_dstrain_var", "Material property name with derivative of G0_pos with strain");
  params.addRequiredParam<MaterialPropertyName>("Emod", "Material property name with Young's Modulus");
  params.addRequiredParam<MaterialPropertyName>("sigmac", "Material property name with strength");
  params.addRequiredCoupledVar("beta", "Auxiliary variable");
  params.addCoupledVar("disp_x", "The x displacement");
  params.addCoupledVar("disp_y", "The y displacement");
  params.addCoupledVar("disp_z", "The z displacement");

  return params;
}

CohesivePFFracBulkRateHistory::CohesivePFFracBulkRateHistory(const InputParameters & parameters):
  Kernel(parameters),
  _gc_prop(getMaterialProperty<Real>("gc_prop_var")),
  _G0_pos(getMaterialProperty<Real>("G0_var")),
 _dG0_pos_dstrain(isParamValid("dG0_dstrain_var") ? &getMaterialProperty<RankTwoTensor>("dG0_dstrain_var"): NULL),
  _Emod(getMaterialProperty<Real>("Emod")),
  _sigmac(getMaterialProperty<Real>("sigmac")),
  _u_old(valueOld()),
  _grad_u_old(gradientOld()),
  _xdisp_coupled(isCoupled("disp_x")),
  _ydisp_coupled(isCoupled("disp_y")),
  _zdisp_coupled(isCoupled("disp_z")),
  _xdisp_var(_xdisp_coupled ? coupled("disp_x") : 0),
  _ydisp_var(_ydisp_coupled ? coupled("disp_y") : 0),
  _zdisp_var(_zdisp_coupled ? coupled("disp_z") : 0),
  _ifOld(getParam<bool>("ifOld")),
  _l(getParam<Real>("l")),
  _p(getParam<Real>("p")),
 _visco(getParam<Real>("visco"))
{
}

Real
CohesivePFFracBulkRateHistory::computeDFDOP(PFFunctionType type)
{
      Real gc = _gc_prop[_qp];

      Real _damage, _beta;

      if (_ifOld){
          _damage = _u_old[_qp];
          _beta   = _grad_test[_i][_qp] * _grad_u_old[_qp];
      }else{
          _damage = _u[_qp];
          _beta   = _grad_test[_i][_qp] * _grad_u[_qp];
      }


      Real _c = 0.375 * _l * gc;
      Real _k = 0.75  * gc / _l;
      Real _m = 1.5 * _Emod[_qp] * gc / (_sigmac[_qp]*_sigmac[_qp]*_l);
      Real _a = (1.0-_damage)*(1.0-_damage);
      Real _b = 1.0 + (_m-2.0)*_damage + (1.0+_p*_m)*_damage*_damage;
      Real _da_dphi = - 2.0 * ( 1.0 - _damage );
      Real _db_dphi = (_m-2.0) + (1.0+_p*_m)*2.0*_damage;

      Real _dg_dphi = _da_dphi/_b - _a/(_b*_b) * _db_dphi;

      Real _psi_e = std::max(_G0_pos[_qp],_k/_m);

      Real x =  - _c * _beta - _test[_i][_qp] * (_k + _dg_dphi * _psi_e );


  switch (type)
  {
    case Residual:
    {
      //if ( x <= 0.0 ) return 0.0;
      return -x/_visco;

    }
    case Jacobian:
    {
      if (_ifOld ) return 0.0;

    	Real _da_dphi2 = 2.0;
      Real _db_dphi2 = 2.0*(1.0+_p*_m);

     	Real _dg_dphi2 = (_da_dphi2*_db_dphi - _db_dphi2*_da_dphi)/(_b*_b) - 2.0*_db_dphi/_b*_dg_dphi;

    	return ( _c * _grad_test[_i][_qp] * _grad_phi[_j][_qp] + _test[_i][_qp] * _dg_dphi2 * _psi_e * _phi[_j][_qp] )/_visco;
    }
    default:
      mooseError("PFFracBulkRate: Invalid type passed - case must be either Residual or Jacobian");
  }
  mooseError("PFFracBulkRate: Invalid type passed");
}

Real
CohesivePFFracBulkRateHistory::computeQpResidual()
{
  return computeDFDOP(Residual);
}

Real
CohesivePFFracBulkRateHistory::computeQpJacobian()
{
  return computeDFDOP(Jacobian);
}



Real
CohesivePFFracBulkRateHistory::computeQpOffDiagJacobian(unsigned int jvar)
{

  if (_ifOld)
      return 0.0;

  unsigned int c_comp;
  bool disp_flag = false;

 // Real x = _l * _betaval[_qp] + 2.0*(1.0-c) * (_G0_pos[_qp]/gc) - c/_l;
 // Real signx = x > 0.0 ? 1.0 : -1.0;

 // Real xfacbeta = -((signx + 1.0)/2.0) / _visco * _l;
 // Real xfac = -((signx + 1.0)/2.0) / _visco * 2.0 * (1.0 - c) / gc;

  Real gc = _gc_prop[_qp];

  Real _damage = _u[_qp];
  Real _c = 0.375 * _l * gc;
  Real _k = 0.75  * gc / _l;

  Real _m = 1.5 * _Emod[_qp] * gc / (_sigmac[_qp]*_sigmac[_qp]*_l);
  Real _a = (1.0-_damage)*(1.0-_damage);
  Real _b = 1.0 + (_m-2.0)*_damage + (1.0+_p*_m)*_damage*_damage;
  Real _da_dphi = - 2.0 * ( 1.0 - _damage );
  Real _db_dphi = (_m-2.0) + (1.0+_p*_m)*2.0*_damage;

  Real _dg_dphi = _da_dphi/_b - _a/(_b*_b) * _db_dphi;
  Real _psi_e = std::max(_G0_pos[_qp],_k/_m);

  //Real x = _c * _betaval[_qp] - _k - _dg_dphi * _psi_e;

  Real x =  - _c * _grad_test[_i][_qp] * _grad_u[_qp] - _test[_i][_qp] * (_k + _dg_dphi * _psi_e );

  Real xfac = 0.0;

  if (_G0_pos[_qp] > _k/_m)
  	xfac = _dg_dphi/_visco;


  //if (jvar == _beta_var)
    //Contribution of auxiliary variable to off diag Jacobian of c
  //  return xfacbeta * _phi[_j][_qp] * _test[_i][_qp];
  if (_xdisp_coupled && jvar == _xdisp_var)
  {
    c_comp = 0;
    disp_flag = true;
  }
  else if (_ydisp_coupled && jvar == _ydisp_var)
  {
    c_comp = 1;
    disp_flag = true;
  }
  else if (_zdisp_coupled && jvar == _zdisp_var)
  {
    c_comp = 2;
    disp_flag = true;
  }
  else
    return 0.0;
  //Contribution of displacements to off diag Jacobian of c
  if (disp_flag && _dG0_pos_dstrain != NULL)
  {
    Real val = 0.0;

    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      val +=  ((*_dG0_pos_dstrain)[_qp](c_comp,i) + (*_dG0_pos_dstrain)[_qp](i,c_comp))/2.0 * _grad_phi[_j][_qp](i);

    return xfac * val * _test[_i][_qp];
  }

  return 0.0;
}
