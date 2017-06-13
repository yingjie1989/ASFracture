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
#ifndef ASFRACTURE_H
#define ASFRACTURE_H

#include "MooseApp.h"

class ASFracture;

template<>
InputParameters validParams<ASFracture>();

class ASFracture : public MooseApp
{
public:
  ASFracture(InputParameters parameters);
  virtual ~ASFracture();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax& syntax, ActionFactory & action_factory);
};

#endif /* ASFRACTURE_H */
