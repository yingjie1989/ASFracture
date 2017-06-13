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

#ifndef MONOPOLEDIRAC_H
#define MONOPOLEDIRAC_H

#include "DiracKernel.h"

class MonopoleDirac;

template<>
InputParameters validParams<MonopoleDirac>();


class MonopoleDirac : public DiracKernel
{
public:
  MonopoleDirac(const InputParameters & parameters);
  
  virtual void addPoints() override;
  virtual Real computeQpResidual() override;

protected:

  Point _point;
  int _dim;
  Real _upcoeff;
  Real _downcoeff;
  Real _rho;
  Real _fL;
  Real _tRT;
  Real _t1;
  Real _tL;
  Real _tP;
  Real _p0;
  Real _d1;

};


#endif /* SOURCEMONOPOLE_H */
