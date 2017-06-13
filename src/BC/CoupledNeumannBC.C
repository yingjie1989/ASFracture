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

#include "CoupledNeumannBC.h"

template<>
InputParameters validParams<CoupledNeumannBC>()
{
	InputParameters params = validParams<IntegratedBC>();
        params.addParam<Real>("Reyold", 1.0, "Value multiplied by the coupled value on the boundary");
        params.addParam<Real>("tol",0.0, "indicate tolerance");
        params.addCoupledVar("some_var",0, "Flux Value at the Boundary");

        return params;
}

CoupledNeumannBC::CoupledNeumannBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _Re(getParam<Real>("Reyold")),
    _tol(getParam<Real>("tol")),
   _some_var_val(coupledValue("some_var"))

{}

Real
CoupledNeumannBC::computeQpResidual()
{
	  return -_test[_i][_qp]*_Re*_some_var_val[_qp];
}
