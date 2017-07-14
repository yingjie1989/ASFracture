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
#ifndef WEIBULLBULKRATEMATERIAL_H
#define WEIBULLBULKRATEMATERIAL_H

#include "Material.h"


//Forward Declarations
class WeibullBulkRateMaterial;

template<>
InputParameters validParams<WeibullBulkRateMaterial>();

/**
 * Empty material for use in simple applications that don't need material properties.
 */
class WeibullBulkRateMaterial : public Material
{
public:
  WeibullBulkRateMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  MaterialProperty<Real> & _weibull;
  MaterialProperty<Real> & _weibull_old;

private:
  int _weibull_modulus;
  Real _specimen_volume;
  Real _Em;
  Real _L0;
  Real _specimen_material_property;
  Real _eta;
  Real _gc;
};

#endif //WEIBULLMATERIAL_H
