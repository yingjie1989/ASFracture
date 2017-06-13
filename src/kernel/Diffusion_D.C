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

#include "Diffusion_D.h"


template<>
InputParameters validParams<Diffusion_D>()
{
  InputParameters params = validParams<Diffusion>();
  params.addParam<MaterialPropertyName>("Diffusivity", "D", "The diffusivity used with the kernel");
  params.addRequiredCoupledVar("fieldvar", "Auxiliary field variable");
  params.addParam<bool>("expdyn",false,"indicate if we are solving explicit dynamics problem");

 return params;
}


Diffusion_D::Diffusion_D(const InputParameters & parameters) :
   Diffusion(parameters),
   _grad_f(coupledGradient("fieldvar")),
   _grad_f_old(coupledGradientOld("fieldvar")),
   _ifOld(getParam<bool>("expdyn")),
   _Diffusivity(getMaterialProperty<Real>("Diffusivity"))


{}

Real
Diffusion_D::computeQpResidual()
{
 if(_ifOld)
	  return _Diffusivity[_qp] * _grad_f_old[_qp] * _grad_test[_i][_qp];

  return _Diffusivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
Diffusion_D::computeQpJacobian()
{
 if(_ifOld)
 	return 0.0;

  return _Diffusivity[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}



