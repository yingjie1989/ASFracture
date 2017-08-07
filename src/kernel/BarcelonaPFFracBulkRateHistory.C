/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "BarcelonaPFFracBulkRateHistory.h"

template<>
InputParameters validParams<BarcelonaPFFracBulkRateHistory>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Kernel to compute bulk energy contribution to damage order parameter residual equation");
  params.addRequiredParam<Real>("l","Interface width");
   params.addRequiredParam<Real>("p","p parameter which influences the cohesive traction separation law");
  params.addRequiredParam<MaterialPropertyName>("gc_prop_var", "Material property name with gc value");
  params.addRequiredParam<MaterialPropertyName>("G0_var", "Material property name with undamaged strain energy driving damage (G0_pos)");
  params.addParam<MaterialPropertyName>("dG0_dstrain_var", "Material property name with derivative of G0_pos with strain");
  params.addRequiredParam<MaterialPropertyName>("Emod", "Material property name with Young's Modulus");
  params.addRequiredParam<MaterialPropertyName>("sigmac", "Material property name with strength");
  params.addCoupledVar("disp_x", "The x displacement");
  params.addCoupledVar("disp_y", "The y displacement");
  params.addCoupledVar("disp_z", "The z displacement");

  return params;
}

BarcelonaPFFracBulkRateHistory::BarcelonaPFFracBulkRateHistory(const InputParameters & parameters):
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
  _l(getParam<Real>("l")),
  _p(getParam<Real>("p"))
{
}

Real
BarcelonaPFFracBulkRateHistory::computeDFDOP(PFFunctionType type)
{
      Real gc = _gc_prop[_qp];

      Real _damage, _beta;

          _damage = _u[_qp];
          _beta   = _grad_test[_i][_qp] * _grad_u[_qp];
          
      Real _c = 0.375 * _l * gc;
      Real _k = 0.75  * gc / _l;
      Real _m = 1.5 * _Emod[_qp] * gc / (_sigmac[_qp]*_sigmac[_qp]*_l);
      Real _a = (1.0-_damage)*(1.0-_damage);
      Real _b = 1.0 + (_m-2.0)*_damage + (1.0+_p*_m)*_damage*_damage;
      Real _da_dphi = - 2.0 * ( 1.0 - _damage );
      Real _db_dphi = (_m-2.0) + (1.0+_p*_m)*2.0*_damage;

      Real _dg_dphi = _da_dphi/_b - _a/(_b*_b) * _db_dphi;

      Real _psi_e = _G0_pos[_qp];

      Real x = 0.5*_l*_l * _beta + _test[_i][_qp] * ( _dg_dphi*_psi_e*_l/gc*4.0/3.0 + 1 );


  switch (type)
  {
    case Residual:
    {
      //if ( x <= 0.0 ) return 0.0;
      return x/_visco;

    }
    case Jacobian:
    {
    	Real _da_dphi2 = 2.0;
      Real _db_dphi2 = 2.0*(1.0+_p*_m);

     	Real _dg_dphi2 = (_da_dphi2*_db_dphi - _db_dphi2*_da_dphi)/(_b*_b) - 2.0*_db_dphi/_b*_dg_dphi;

    	return ( 0.5*_l*_l*_grad_test[_i][_qp]*_grad_phi[_j][_qp] + _test[_i][_qp] * _dg_dphi2 * _phi[_j][_qp]*_psi_e*_l/gc*4.0/3.0 );
    }
    default:
      mooseError("PFFracBulkRate: Invalid type passed - case must be either Residual or Jacobian");
  }
  mooseError("PFFracBulkRate: Invalid type passed");
}

Real
BarcelonaPFFracBulkRateHistory::computeQpResidual()
{
  return computeDFDOP(Residual);
}

Real
BarcelonaPFFracBulkRateHistory::computeQpJacobian()
{
  return computeDFDOP(Jacobian);
}



Real
BarcelonaPFFracBulkRateHistory::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
