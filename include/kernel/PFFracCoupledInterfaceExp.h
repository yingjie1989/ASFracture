/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PFFRACCOUPLEDINTERFACEEXP_H
#define PFFRACCOUPLEDINTERFACEEXP_H

#include "KernelGrad.h"

/**
 * Phase-field fracture model
 * This class computes the contribution to residual and jacobian of the variable beta
 * by the grad of c (damage variable)
 * Refer to formulation: Miehe et. al., Int. J. Num. Methods Engg., 2010, 83. 1273-1311 Equation 63
 */

//Forward Declarations
class PFFracCoupledInterfaceExp;

template<>
InputParameters validParams<PFFracCoupledInterfaceExp>();

class PFFracCoupledInterfaceExp : public KernelGrad
{
public:

  PFFracCoupledInterfaceExp(const InputParameters & parameters);

protected:
  virtual RealGradient precomputeQpResidual();
  virtual RealGradient precomputeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableValue & _c;
  const VariableGradient & _grad_c;
  unsigned int _c_var;

private:
};
#endif //PFFRACCOUPLEDINTERFACEEXP_H
