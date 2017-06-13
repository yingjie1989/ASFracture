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

#ifndef DIFFUSION_D_H
#define DIFFUSION_D_H

#include "Diffusion.h"

class Diffusion_D;

template<>
InputParameters validParams<Diffusion_D>();


class Diffusion_D : public Diffusion
{
public:
  Diffusion_D(const InputParameters & parameters);
  const VariableGradient & _grad_f;
  const VariableGradient & _grad_f_old;
  bool _ifOld;
  const MaterialProperty<Real> & _Diffusivity;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
};


#endif /* DIFFUSION_D_H */
