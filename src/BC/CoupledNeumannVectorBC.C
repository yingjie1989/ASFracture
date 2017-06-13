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

#include "CoupledNeumannVectorBC.h"

template<>
InputParameters validParams<CoupledNeumannVectorBC>()
{
	InputParameters params = validParams<IntegratedBC>();
        params.addParam<Real>("Reyold", 1.0, "Value multiplied by the coupled value on the boundary");
        params.addParam<Real>("tol",0.0, "indicate tolerance");
        params.addCoupledVar("some_var_x",0, "Flux_x Value at the Boundary");
	params.addCoupledVar("some_var_y",0, "Flux_y Value at the Boundary");
	params.addCoupledVar("some_var_z",0, "Flux_z Value at the Boundary");


        return params;
}

CoupledNeumannVectorBC::CoupledNeumannVectorBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _Re(getParam<Real>("Reyold")),
    _tol(getParam<Real>("tol")),
    _some_var_x(coupledValue("some_var_x")),
    _some_var_y(coupledValue("some_var_y")),
    _some_var_z(coupledValue("some_var_z"))

{}

Real
CoupledNeumannVectorBC::computeQpResidual()
{
	  return -_test[_i][_qp] * _Re * ( _some_var_x[_qp] * _normals[_qp](0) + _some_var_y[_qp] * _normals[_qp](1) + _some_var_z[_qp] * _normals[_qp](2) ) ;
}
