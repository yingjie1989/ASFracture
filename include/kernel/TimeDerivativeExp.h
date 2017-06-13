/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef TIMEDERIVATIVEEXP_H
#define TIMEDERIVATIVEEXP_H

#include "Kernel.h"
#include "Material.h"

//Forward Declarations
class TimeDerivativeExp;

template<>
InputParameters validParams<TimeDerivativeExp>();

class TimeDerivativeExp : public Kernel
{
public:

  TimeDerivativeExp(const InputParameters & parameters);

  virtual void computeJacobian() override;


protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

private:
  Real _coeff;
  bool _lumped;
  const VariableValue & _u_old;
  const VariableValue & _u_nodal;
  const VariableValue & _u_nodal_old;

  };

#endif //TIMEDERIVATIVEEXP_H
