/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ElasticEnergyMaterial2DFrac_H
#define ElasticEnergyMaterial2DFrac_H

#include "DerivativeFunctionMaterialBase.h"

// Forward Declaration
class ElasticEnergyMaterial2DFrac;
class RankTwoTensor;
class RankFourTensor;

template<>
InputParameters validParams<DerivativeFunctionMaterialBase>();

/**
 * Material class to compute the elastic free energy and its derivatives
 */
class ElasticEnergyMaterial2DFrac : public DerivativeFunctionMaterialBase
{
public:
  ElasticEnergyMaterial2DFrac(const InputParameters & parameters);

protected:
  virtual Real computeF();
  virtual Real computeDF(unsigned int);
  virtual Real computeD2F(unsigned int, unsigned int);

  std::string _base_name;

  const MaterialProperty<Real> & _energy_tensile;
  const MaterialProperty<Real> & _d_energy_tensile_dc;
  const MaterialProperty<Real> & _d2_energy_tensile_dc2;


};

#endif //ElasticEnergyMaterial2DFrac_H
