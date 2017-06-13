/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PFFracRandomBulkRateMaterial.h"

template<>
InputParameters validParams<PFFracRandomBulkRateMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Material properties used in phase-field fracture damage evolution kernel");
  params.addParam<FunctionName>("function", "Function describing energy release rate type parameter distribution");
  params.addCoupledVar("beta",0.0, "perturbation variable");
  params.addParam<Real>("gc", 1.0, "Energy release rate type parameter");
  params.addParam<Real>("pC", 0.0, "Perturbation of Energy release rate");
  
  return params;
}

PFFracRandomBulkRateMaterial::PFFracRandomBulkRateMaterial(const InputParameters & parameters) :
    Material(parameters),
    _gc(getParam<Real>("gc")),
    _perturbCoeff(getParam<Real>("pC")),
    _betaval(coupledValue("beta")), 
    _gc_prop(declareProperty<Real>("gc_prop")),
    _gc_prop_old(declarePropertyOld<Real>("gc_prop")),
    _function_prop(isParamValid("function") ? &getFunction("function") : NULL)
{

    setRandomResetFrequency(EXEC_INITIAL); 
}

void
PFFracRandomBulkRateMaterial::initQpStatefulProperties()
{
     Real random_real = getRandomReal();

      _gc_prop[_qp] = (1.0 + _betaval[_qp] - _perturbCoeff + 2*_perturbCoeff*random_real) * _gc;
   
      _gc_prop_old[_qp] = _gc_prop[_qp]; 
}

void
PFFracRandomBulkRateMaterial::computeQpProperties()
{
}

void
PFFracRandomBulkRateMaterial::getProp()
{
  if (_function_prop != NULL)
    _gc_prop[_qp] = _function_prop->value(_t, _q_point[_qp]);
}

