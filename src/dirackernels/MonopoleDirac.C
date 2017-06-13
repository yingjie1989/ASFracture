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

#include "MonopoleDirac.h"
# define _PI 3.14159265358979323846  /* pi */

template<>
InputParameters validParams<MonopoleDirac>()
{
  InputParameters params = validParams<DiracKernel>();
  params.addRequiredParam<Point>("point", "The x,y,z coordinates of the point"); 
  params.addRequiredParam<int>("dim", "dimension of problem");
  params.addRequiredParam<Real>("upcoeff", "upcoefficient");
  params.addRequiredParam<Real>("downcoeff", "downcoefficient");
  params.addRequiredParam<Real>("rho","density");
  params.addRequiredParam<Real>("fL","fLcoefficient");
  params.addRequiredParam<Real>("tRT","tRTcoefficient");
  params.addRequiredParam<Real>("t1","t1coefficient");
  params.addRequiredParam<Real>("tL","tLcoefficient");
  params.addRequiredParam<Real>("tP","tPcoefficient");
  params.addRequiredParam<Real>("p0","p0coefficient");
  params.addRequiredParam<Real>("d1","d1coefficient");
  return params;
}


MonopoleDirac::MonopoleDirac(const InputParameters & parameters) :
   DiracKernel(parameters),
    _point(getParam<Point>("point")),
    _dim(getParam<int>("dim")), 
    _upcoeff(getParam<Real>("upcoeff")),
    _downcoeff(getParam<Real>("downcoeff")),
    _rho(getParam<Real>("rho")),
    _fL(getParam<Real>("fL")),
    _tRT(getParam<Real>("tRT")),
    _t1(getParam<Real>("t1")),
    _tL(getParam<Real>("tL")),
    _tP(getParam<Real>("tP")),
    _p0(getParam<Real>("p0")),
    _d1(getParam<Real>("d1"))
{}

void
MonopoleDirac::addPoints()
{
  // Add a point from the input file
     addPoint(_point);
}

Real
MonopoleDirac::computeQpResidual()
{ 

    //Real _rad = std::sqrt( (_q_point[_qp](0) - _point(0))*(_q_point[_qp](0) - _point(0)) + (_q_point[_qp](1) - _point(1))*(_q_point[_qp](1) - _point(1)) + (_q_point[_qp](2) - _point(2))*(_q_point[_qp](2) - _point(2)) );

    Real PseudoTime = _t;

    while(PseudoTime > _tP){

	PseudoTime += - _tP;
    }

    Real func1 = _upcoeff/_downcoeff * (1+std::tanh((PseudoTime - _t1)/_tRT ))*std::exp(-(PseudoTime - _t1)/_tL)*std::cos(2*_PI*_fL*(PseudoTime - _t1) + _PI/3.0);

    Real _value = 0.0;
      
    _value = 2 * (_dim - 1) *_PI/_rho * std::max(func1,0.0) * _p0 * _d1; 

    return -_test[_i][_qp] * _value;
}







