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

#ifndef COUPLEDNEUMANNVECTORBC_H
#define COUPLEDNEUMANNVECTORBC_H

#include "IntegratedBC.h"

//Forward Declarations

class CoupledNeumannVectorBC;

template<>
InputParameters validParams<CoupledNeumannVectorBC>();

class CoupledNeumannVectorBC : public IntegratedBC
{
public:


 CoupledNeumannVectorBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:


 Real _Re;
  Real _tol;

 const VariableValue & _some_var_x;
 const VariableValue & _some_var_y;
 const VariableValue & _some_var_z;


};

#endif //COUPLEDNEUMANNBCVECTOR_H

