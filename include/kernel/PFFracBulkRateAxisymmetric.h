/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PFFRACBULKRATEAXISYMMETRIC_H
#define PFFRACBULKRATEAXISYMMETRIC_H
/**
 * Phase field based fracture model
 * This kernel computes the residual and jacobian for bulk free energy contribution to cU
 * Refer to Formulation: Miehe et. al., Int. J. Num. Methods Engg., 2010, 83. 1273-1311 Equation 63
 */
#include "Kernel.h"
#include "RankTwoTensor.h"

//Forward Declarations
class PFFracBulkRateAxisymmetric;

template<>
InputParameters validParams<PFFracBulkRateAxisymmetric>();

class PFFracBulkRateAxisymmetric : public Kernel
{
public:

  PFFracBulkRateAxisymmetric(const InputParameters & parameters);

protected:


  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  ///Critical energy release rate for fracture
  const MaterialProperty<Real> & _gc_prop;
  ///Contribution of umdamaged strain energy to damage evolution
  const MaterialProperty<Real> & _G0_pos;
  ///Variation of undamaged strain energy driving damage evolution with strain
  const MaterialProperty<RankTwoTensor> * _dG0_pos_dstrain;
  ///Auxiliary variable: beta = Laplacian of c
  const VariableValue & _betaval;
  const unsigned int _beta_var;

  const bool _xdisp_coupled;
  const bool _ydisp_coupled;
  const bool _zdisp_coupled;

  const unsigned int _xdisp_var;
  const unsigned int _ydisp_var;
  const unsigned int _zdisp_var;

  ///Characteristic length, controls damage zone thickness
  Real _l;
  Real _kdamage;
  ///Viscosity parameter ( visco -> 0, rate independent )
  Real _visco;

 private:

};
#endif //PFFRACBULKRATEAXISYMMETRIC_H
