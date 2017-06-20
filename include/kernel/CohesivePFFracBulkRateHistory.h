/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COHESIVEPFRACBULKRATEHSITORY_H
#define COHESIVEPFRACBULKRATEHSITORY_H
/**
 * Phase field based fracture model
 * This kernel computes the residual and jacobian for bulk free energy contribution to c
 * Refer to Formulation: Miehe et. al., Int. J. Num. Methods Engg., 2010, 83. 1273-1311 Equation 63
 */
#include "Kernel.h"
#include "RankTwoTensor.h"

//Forward Declarations
class CohesivePFFracBulkRateHistory;

template<>
InputParameters validParams<CohesivePFFracBulkRateHistory>();

class CohesivePFFracBulkRateHistory : public Kernel
{
public:

  CohesivePFFracBulkRateHistory(const InputParameters & parameters);

protected:

  enum PFFunctionType
  {
    Jacobian,
    Residual
  };

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  virtual Real computeDFDOP(PFFunctionType type);
  ///Critical energy release rate for fracture
  const MaterialProperty<Real> & _gc_prop;
  ///Contribution of umdamaged strain energy to damage evolution
  const MaterialProperty<Real> & _G0_pos;
  ///Variation of undamaged strain energy driving damage evolution with strain
  const MaterialProperty<RankTwoTensor> * _dG0_pos_dstrain;
  //Young's Modulus
  const MaterialProperty<Real> & _Emod;
  //Strenth
  const MaterialProperty<Real> & _sigmac;
  ///Auxiliary variable: beta = Laplacian of c
  const VariableValue & _u_old;
  const VariableGradient & _grad_u_old;
  const bool _xdisp_coupled;
  const bool _ydisp_coupled;
  const bool _zdisp_coupled;

  const unsigned int _xdisp_var;
  const unsigned int _ydisp_var;
  const unsigned int _zdisp_var;

  bool _ifOld;
  ///Characteristic length, controls damage zone thickness
  Real _l;
  //p parameter
  Real _p;

 ///Viscosity parameter ( visco -> 0, rate independent )
  Real _visco;

 private:

};
#endif //CohesivePFFracBulkRateHistory_H
