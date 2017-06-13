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

#ifndef COUPLEDNEUMANNBC_H
#define COUPLEDNEUMANNBC_H

#include "IntegratedBC.h"

//Forward Declarations

class CoupledNeumannBC;

template<>
InputParameters validParams<CoupledNeumannBC>();

class CoupledNeumannBC : public IntegratedBC
{
public:


 CoupledNeumannBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:


 Real _Re;
  Real _tol;

 const VariableValue & _some_var_val;

};

#endif //COUPLEDNEUMANNBC_H

