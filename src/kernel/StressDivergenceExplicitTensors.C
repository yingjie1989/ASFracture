/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergenceExplicitTensors.h"
#include "Assembly.h"
#include "ElasticityTensorTools.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"


template<>
InputParameters validParams<StressDivergenceExplicitTensors>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Stress divergence kernel for phase-field fracture: Additionally computes off diagonal damage dependent Jacobian components");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addRequiredParam<unsigned int>("component",
                                            "An integer corresponding to the direction "
                                            "the variable this kernel acts in. (0 for x, "
                                            "1 for y, 2 for z)");
  //params.addRequiredParam<MaterialPropertyName>("gc_prop_var", "Material property name with gc value");

  return params;
}


StressDivergenceExplicitTensors::StressDivergenceExplicitTensors(const InputParameters & parameters) :
    Kernel(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _grad_disp(3),

    //_strain_old(declareProperty<RankTwoTensor>("mechanical_strain")),
    //_stress_old(declareProperty<RankTwoTensor>("stress_old")),
    _elasticity_tensor_name("elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _component(getParam<unsigned int>("component"))

{
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp[i] = &coupledValueOld("displacements", i);
    _grad_disp[i] = &coupledGradientOld("displacements", i);
  }

  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _disp[i] = &_zero;
    _grad_disp[i] = &_grad_zero;
  }

}

/*void
StressDivergenceExplicitTensors::computeStrain()
{

}

void
StressDivergenceExplicitTensors::computeStress()
{
}*/

Real
StressDivergenceExplicitTensors::computeQpResidual()
{

 RankTwoTensor _strain_old;
 RankTwoTensor _stress_old;

 //for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
 //{
   // strain = (grad_disp + grad_disp^T)/2
   RankTwoTensor grad_tensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);

  _strain_old = (grad_tensor + grad_tensor.transpose()) / 2.0;
 //}
  _stress_old = _elasticity_tensor[_qp] * _strain_old;

  Real residual = _stress_old.row(_component) * _grad_test[_i][_qp];

  return residual;
}

Real
StressDivergenceExplicitTensors::computeQpJacobian()
{
	return 0.0;
}


Real
StressDivergenceExplicitTensors::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.0;
}
