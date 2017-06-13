/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "TimeDerivativeExp.h"
#include "SubProblem.h"
#include "Assembly.h"
#include "MooseVariable.h"
// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<TimeDerivativeExp>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the residual for the interial force (M*accel) and the contribution of mass dependent Rayleigh damping and HHT time integration scheme [eta*M*((1+alpha)vel-alpha*vel_old)]");
  params.set<bool>("use_displaced_mesh") = true;
  params.addParam<Real>("coeff",1.0,"coefficient in front of time derivative");
  params.addParam<bool>("use_lumped_mass",false,"indicate whether use lumped mass matrix");
  return params;
}

TimeDerivativeExp::TimeDerivativeExp(const InputParameters & parameters) :
    Kernel(parameters),
    _coeff(getParam<Real>("coeff")),
    _lumped(getParam<bool>("use_lumped_mass")),
    _u_old(valueOld()),
    _u_nodal(_var.nodalValue()),
    _u_nodal_old(_var.nodalValueOld())
{}

Real
TimeDerivativeExp::computeQpResidual()
{

  Real vel = 0.0;
  if (_lumped)
    vel += 1./_dt * ( _u_nodal[_qp] - _u_nodal_old[_qp] );
  else
    vel += 1./_dt * ( _u[_qp] - _u_old[_qp] );

  return _test[_i][_qp] * _coeff * vel;

}

Real
TimeDerivativeExp::computeQpJacobian()
{
  if (_lumped)
    return _test[_i][_qp] * _coeff / _dt;
 else
    return _test[_i][_qp] * _coeff / _dt * _phi[_j][_qp];
}

void
TimeDerivativeExp::computeJacobian()
{
  if (_lumped){
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  for (_i = 0; _i < _test.size(); _i++)
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
    }else{

      Kernel::computeJacobian();
    }
}
