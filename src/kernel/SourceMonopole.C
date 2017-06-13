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

#include "SourceMonopole.h"
# define _PI 3.14159265358979323846  /* pi */

template<>
InputParameters validParams<SourceMonopole>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<std::vector<Real> >("coord", "position of source");
  params.addRequiredParam<Real>("size", "size of source");
  params.addRequiredParam<Real>("upcoeff", "size of source");
  params.addRequiredParam<Real>("downcoeff", "size of source");
  params.addRequiredParam<Real>("fL","fLcoefficient");
  params.addRequiredParam<Real>("tRT","tRTcoefficient");
  params.addRequiredParam<Real>("t1","t1coefficient");
  params.addRequiredParam<Real>("tL","tLcoefficient");
  params.addRequiredParam<Real>("p0","p0coefficient");
  params.addRequiredParam<Real>("d1","d1coefficient");
  params.addRequiredParam<Real>("rho_c","d1coefficient");
  return params;
}


SourceMonopole::SourceMonopole(const InputParameters & parameters) :
   Kernel(parameters),
    _coord(getParam<std::vector<Real> >("coord")),
    _size(getParam<Real>("size")),
    _upcoeff(getParam<Real>("upcoeff")),
    _downcoeff(getParam<Real>("downcoeff")),
    _fL(getParam<Real>("fL")),
    _tRT(getParam<Real>("tRT")),
    _t1(getParam<Real>("t1")),
    _tL(getParam<Real>("tL")),
    _p0(getParam<Real>("p0")),
    _d1(getParam<Real>("d1")),
    _rho_c(getParam<Real>("rho_c"))
{}

Real
SourceMonopole::computeQpResidual()
{ 
    Real func1 = 0.0;

    Real _rad = std::sqrt( (_q_point[_qp](0) - _coord[0])*(_q_point[_qp](0) - _coord[0]) + (_q_point[_qp](1) - _coord[1])*(_q_point[_qp](1) - _coord[1]) + (_q_point[_qp](2) - _coord[2])*(_q_point[_qp](2) - _coord[2]) );

    if(_rad <= _size)
 	func1 += _upcoeff/_downcoeff * (1+std::tanh((_t - _t1)/_tRT ))*std::exp(-(_t - _t1)/_tL)*std::cos(2*_PI*_fL*(_t-_t1) + _PI/3.0);

    Real func = 4*_PI/_rho_c * std::max(func1,0.0) * _p0 * _d1; 
	
    return -_test[_i][_qp] * func;
}

Real
SourceMonopole::computeQpJacobian()
{
	return 0.0;
}







