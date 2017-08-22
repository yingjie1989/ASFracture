/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef LINEARELASTICPFDAMAGE_H
#define LINEARELASTICPFDAMAGE_H

#include "ComputeStressBase.h"
#include "Function.h"
#include "RankTwoTensor.h"

/**
 * Phase-field fracture
 * This class computes the energy contribution to damage growth
 * Small strain Isotropic Elastic formulation
 * Stiffness matrix scaled for heterogeneous elasticity property
 */
class LinearElasticPFDamage : public ComputeStressBase
{
public:
  LinearElasticPFDamage(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;
//  virtual void updateVar();
//  virtual void updateJacobian();

  bool _ifenergy;
  bool _ifstress;
  bool _ifstrain;
  const VariableValue & _c;
  /// Small number to avoid non-positive definiteness at or near complete damage
  Real _kdamage;

  MaterialProperty<Real> & _G0_pos;
  MaterialProperty<Real> & _G0_pos_old;
  MaterialProperty<RankTwoTensor> & _dstress_dc;
  MaterialProperty<RankTwoTensor> & _dG0_pos_dstrain;

  std::vector<RankTwoTensor> _etens;
  std::vector<Real> _epos;
  std::vector<Real> _eigval;
  RankTwoTensor _eigvec;
};

#endif //LINEARELASTICPFDAMAGE_H
