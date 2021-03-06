/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STRESSDIVERGENCERZPFFRACTENSORS_H
#define STRESSDIVERGENCERZPFFRACTENSORS_H

#include "StressDivergenceRZTensors.h"
#include "Material.h"
#include "DerivativeMaterialInterface.h"

/**
 * This class computes the off-diagonal Jacobian component of stress divergence residual system
 * Contribution from damage order parameter c
 * Residual calculated in StressDivergenceTensors
 * Useful if user wants to add the off diagonal Jacobian term
 */

class StressDivergenceRZPFFracTensors;

template<>
InputParameters validParams<StressDivergenceRZPFFracTensors>();

class StressDivergenceRZPFFracTensors : public DerivativeMaterialInterface<StressDivergenceRZTensors>
{
public:
  StressDivergenceRZPFFracTensors(const InputParameters & parameters);

protected:
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const bool _c_coupled;
  const unsigned int _c_var;
  const MaterialProperty<RankTwoTensor> & _d_stress_dc;
};

#endif //STRESSDIVERGENCERZPFFRACTENSORS_H
