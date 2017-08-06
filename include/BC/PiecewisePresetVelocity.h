/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PIECEWISEPRESETVELOCITY_H
#define PIECEWISEPRESETVELOCITY_H

#include "PiecewisePresetNodalBC.h"

class PiecewisePresetVelocity : public PresetNodalBC
{
public:
  PiecewisePresetVelocity(const InputParameters & parameters);

protected:
  virtual Real computeQpValue();

  const VariableValue & _u_old;
  const Real _t0;
  const Real _velocity;
  Function & _function;
};

template <>
InputParameters validParams<PresetVelocity>();

#endif /* PIECEWISEPRESETVELOCITY_H */
