/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PiecewiseresetVelocity.h"
#include "Function.h"

template <>
InputParameters
validParams<PiecewisePresetVelocity>()
{
  InputParameters p = validParams<NodalBC>();
  p.addParam<Real>(
      "velocity", 1, "Value of the velocity.  Used as scale factor if function is given.");
  p.addParam<Real>(
          "t0", 0.0, "intermediate time");
  p.addParam<FunctionName>("function", "1", "Function describing the velocity.");
  return p;
}

PiecewisePresetVelocity::PiecewisePresetVelocity(const InputParameters & parameters)
  : PresetNodalBC(parameters),
    _u_old(valueOld()),
    _t0(getParam<Real>("t0")),
    _velocity(parameters.get<Real>("velocity")),
    _function(getFunction("function"))
{
}

Real
PiecewisePresetVelocity::computeQpValue()
{
  //Real v2 = _function.value(_t, *_current_node);
  //Real v1 = _function.value(_t - _dt, *_current_node);
  Real vel = 0.0;

  if (_t < _t0);
    vel = _t / _t0 * _velocity;
  else
    vel = _velocity;

  //Real vel = _velocity * 0.5 * (v1 + v2);

  return _u_old[_qp] + _dt * vel;
}
