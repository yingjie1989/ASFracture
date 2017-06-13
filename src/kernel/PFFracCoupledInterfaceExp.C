/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PFFracCoupledInterfaceExp.h"

template<>
InputParameters validParams<PFFracCoupledInterfaceExp>()
{
  InputParameters params = validParams<KernelGrad>();
  params.addClassDescription("Phase-field fracture residual for beta variable: Contribution from gradient of damage order parameter");
  params.addRequiredCoupledVar("c", "Order parameter for damage");

  return params;
}

PFFracCoupledInterfaceExp::PFFracCoupledInterfaceExp(const InputParameters & parameters):
  KernelGrad(parameters),
  _c(coupledValueOld("c")),
  _grad_c(coupledGradientOld("c")),
  _c_var(coupled("c"))
{
}

RealGradient
PFFracCoupledInterfaceExp::precomputeQpResidual()
{
  return  _grad_c[_qp];
}

RealGradient
PFFracCoupledInterfaceExp::precomputeQpJacobian()
{
  return 0.0;
}
/**
 * Contributes only to the off-diagonal Jacobian term
 * for variable beta
 */
Real
PFFracCoupledInterfaceExp::computeQpOffDiagJacobian(unsigned int jvar)
{
    return 0.0;
}
