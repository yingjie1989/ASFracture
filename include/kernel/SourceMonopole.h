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

#ifndef SOURCEMONOPOLE_H
#define SOURCEMONOPOLE_H

#include "Kernel.h"

class SourceMonopole;

template<>
InputParameters validParams<Kernel>();


class SourceMonopole : public Kernel
{
public:
  SourceMonopole(const InputParameters & parameters);

  std::vector<Real> _coord;
  Real _size;
  Real _upcoeff;
  Real _downcoeff;
  Real _fL;
  Real _tRT;
  Real _t1;
  Real _tL;
  Real _p0;
  Real _d1;
  Real _rho_c;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
};


#endif /* SOURCEMONOPOLE_H */
