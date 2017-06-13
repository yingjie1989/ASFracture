/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COHESIVELINEARISOELASTICPFDAMAGE_H
#define COHESIVELINEARISOELASTICPFDAMAGE_H

#include "ComputeStressBase.h"
#include "Function.h"

/**
 * Phase-field fracture
 * This class computes the energy contribution to damage growth
 * Small strain Isotropic Elastic formulation
 * Stiffness matrix scaled for heterogeneous elasticity property
 */
class CohesiveLinearIsoElasticPFDamage : public ComputeStressBase
{
public:
  CohesiveLinearIsoElasticPFDamage(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress();
  virtual void updateVar();
  virtual void updateJacobian();

  const VariableValue & _c;
  /// Small number to avoid non-positive definiteness at or near complete damage
  Real _kdamage;
  bool _historyEng;
  const MaterialProperty<Real> & _gc_prop;
  //Young's Modulus
  const MaterialProperty<Real> & _Emod;
  //Strenth
  const MaterialProperty<Real> & _sigmac;
  ///Characteristic length, controls damage zone thickness
  Real _l;
  //p parameter
  Real _p;

  MaterialProperty<Real> & _G0_pos;
  MaterialProperty<Real> & _G0_pos_old;
  MaterialProperty<RankTwoTensor> & _dstress_dc;
  MaterialProperty<RankTwoTensor> & _dG0_pos_dstrain;

  std::vector<RankTwoTensor> _etens;
  std::vector<Real> _epos;
  std::vector<Real> _eigval;
  RankTwoTensor _eigvec;
};

#endif //COHESIVELINEARISOELASTICPFDAMAGE_H
