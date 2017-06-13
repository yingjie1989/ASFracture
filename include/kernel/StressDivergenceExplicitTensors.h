/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STRESSDIVERGENCEEXPLICITTENSORS_H
#define STRESSDIVERGENCEEXPLICITTENSORS_H

#include "StressDivergenceTensors.h"
#include "Material.h"
#include "Kernel.h"

/**
 * This class computes the off-diagonal Jacobian component of stress divergence residual system
 * Contribution from damage order parameter c
 * Residual calculated in StressDivergenceTensors
 * Useful if user wants to add the off diagonal Jacobian term
 */

class StressDivergenceExplicitTensors;

template<>
InputParameters validParams<StressDivergenceExplicitTensors>();

class StressDivergenceExplicitTensors : public Kernel
{
public:
  StressDivergenceExplicitTensors(const InputParameters & parameters);

protected:

  unsigned int _ndisp;
  std::vector<const VariableValue *> _disp;
  std::vector<const VariableGradient *> _grad_disp;

  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  unsigned int _component;

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

};

#endif //STRESSDIVERGENCEEXPLICITTENSORS_H
