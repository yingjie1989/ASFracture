/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PFFRACRANDOMBULKRATEMATERIAL_H
#define PFFRACRANDOMBULKRATEMATERIAL_H

#include "Material.h"
#include "Function.h"

/**
 * Phase-field fracture
 * This class obtains critical energy release rate (gc) value
 * Used by PFFRacBulkRate
 */

class PFFracRandomBulkRateMaterial;

template<>
InputParameters validParams<PFFracRandomBulkRateMaterial>();

class PFFracRandomBulkRateMaterial : public Material
{
public:
  PFFracRandomBulkRateMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();
  /**
   * This function obtains the value of gc
   * Must be overidden by the user for heterogeneous gc
   */
  virtual void getProp();

  ///Input parameter for homogeneous gc
  Real _gc;
  Real _perturbCoeff;
  const VariableValue & _betaval;
  ///Material property where the gc values are stored
  MaterialProperty<Real> &_gc_prop;
  MaterialProperty<Real> &_gc_prop_old;

  ///Function to specify varying gc
  Function * _function_prop;

private:

};

#endif //PFFRACRANDOMBULKRATEMATERIAL_H
