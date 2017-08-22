/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "LinearElasticPFDamage.h"
#include "libmesh/utility.h"

template <> InputParameters validParams<LinearElasticPFDamage>() {
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Phase-field fracture model energy contribution "
                             "to damage growth-isotropic elasticity and "
                             "undamaged stress under compressive strain");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<Real>("kdamage", 1e-6, "Stiffness of damaged matrix");
  params.addParam<bool>("ifenergy", false, "whether to calculate energy");
  params.addParam<bool>("ifstress", false,
                        "whether to calculate stress & jacobian");
  params.addParam<bool>("ifstrain", false,
                        "whether to calculate strain & jacobian");
  return params;
}

LinearElasticPFDamage::LinearElasticPFDamage(const InputParameters &parameters)
    : ComputeStressBase(parameters), _ifenergy(getParam<bool>("ifenergy")),
      _ifstress(getParam<bool>("ifstress")),
      _ifstrain(getParam<bool>("ifstrain")), _c(coupledValue("c")),
      _kdamage(getParam<Real>("kdamage")),
      _G0_pos(declareProperty<Real>("G0_pos")),
      _G0_pos_old(declarePropertyOld<Real>("G0_pos")),
      _dstress_dc(declarePropertyDerivative<RankTwoTensor>(
          _base_name + "stress", getVar("c", 0)->name())),
      _dG0_pos_dstrain(declareProperty<RankTwoTensor>("dG0_pos_dstrain")),
      _etens(LIBMESH_DIM), _epos(LIBMESH_DIM), _eigval(LIBMESH_DIM) {}

void LinearElasticPFDamage::initQpStatefulProperties() {
  ComputeStressBase::initQpStatefulProperties();

  if (_t > 0)
    _G0_pos_old[_qp] = _G0_pos[_qp];
  else {
    _G0_pos[_qp] = 0.0;
    _G0_pos_old[_qp] = 0.0;
  }
  _dG0_pos_dstrain[_qp] = _stress[_qp];
  _dstress_dc[_qp] = -_stress[_qp] * (2.0 * (1.0 - _c[_qp]));
}

void LinearElasticPFDamage::computeQpStress() {
  RankTwoTensor stress0pos, stress0neg, stress0;
  // Isotropic elasticity is assumed
  Real c = _c[_qp];
  Real xfac = (Utility::pow<2>(1.0 - c)) + _kdamage;

  if (_ifstress) {
    stress0 = _elasticity_tensor[_qp] * _mechanical_strain[_qp];

    stress0.symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

    // Tensors of outerproduct of eigen vectors
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
          _etens[i](j, k) = _eigvec(j, i) * _eigvec(k, i);

    for (unsigned int i = 0; i < LIBMESH_DIM; ++i) {
      if (_eigval[i] >= 0.0) {
        stress0pos += _etens[i] * _eigval[i];
      } else
        stress0neg += _etens[i] * _eigval[i];
    }
  }
  if (_ifstrain){
    _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

    // Tensors of outerproduct of eigen vectors
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
          _etens[i](j, k) = _eigvec(j, i) * _eigvec(k, i);

    RankTwoTensor strain0pos, strain0neg;  
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i) {
      if (_eigval[i] >= 0.0) {
        strain0pos += _etens[i] * _eigval[i];
      } else
        strain0neg += _etens[i] * _eigval[i];
    }
    stress0pos = _elasticity_tensor[_qp] * strain0pos;
    stress0neg = _elasticity_tensor[_qp] * strain0neg;
  }
  // Damage associated with positive component of stress
  _stress[_qp] = stress0pos * xfac + stress0neg;

  if (_ifenergy) {

    Real G0_trial = 0.5 * stress0pos.doubleContraction(_mechanical_strain[_qp]);

    if (G0_trial > _G0_pos_old[_qp]) {
      _G0_pos[_qp] = G0_trial;
      _dG0_pos_dstrain[_qp] = stress0pos;
    } else {
      _G0_pos[_qp] = _G0_pos_old[_qp];
      _dG0_pos_dstrain[_qp] = stress0pos * 0.0;
    }
    // Used in StressDivergencePFFracTensors Jacobian
    _dstress_dc[_qp] = -stress0pos * (2.0 * (1.0 - c));
  }

  if (_ifstress || _ifstrain) {
    // Used in StressDivergencePFFracTensors Jacobian
    int alpha, beta, i, j, k, l;
    int dim = LIBMESH_DIM;

    // Compute the second derivatives of the eigenvalues w.r.t. strain tensor
    std::vector<RankFourTensor> d2eigvals(dim);
    RankFourTensor dstress_plus, dstress_minus;
    RankFourTensor _Jacob_plus, _Jacob_minus;

    for (alpha = 0; alpha < dim; ++alpha)
      for (beta = 0; beta < dim; ++beta) {
        if (_eigval[alpha] == _eigval[beta])
          continue;

        for (i = 0; i < dim; ++i)
          for (j = 0; j < dim; ++j)
            for (k = 0; k < dim; ++k)
              for (l = 0; l < dim; ++l) {
                d2eigvals[alpha](i, j, k, l) +=
                    0.5 *
                    (_eigvec(beta, i) * _eigvec(alpha, j) +
                     _eigvec(alpha, i) * _eigvec(beta, j)) *
                    (_eigvec(beta, k) * _eigvec(alpha, l) +
                     _eigvec(beta, l) * _eigvec(alpha, k)) /
                    (_eigval[alpha] - _eigval[beta]);
              }
      }

    for (alpha = 0; alpha < dim; ++alpha) // loop in eigenvalues
    {
      for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
          for (k = 0; k < dim; ++k) {
            for (l = 0; l < dim; ++l) {
              if (_eigval[alpha] >= 0.0)
                dstress_plus(i, j, k, l) +=
                    _eigval[alpha] * d2eigvals[alpha](i, j, k, l) +
                    _etens[alpha](i, j) * _etens[alpha](k, l);
              if (_eigval[alpha] < 0.0)
                dstress_minus(i, j, k, l) +=
                    _eigval[alpha] * d2eigvals[alpha](i, j, k, l) +
                    _etens[alpha](i, j) * _etens[alpha](k, l);
            }
          }
        }
      }
    }

    _Jacobian_mult[_qp] =
        (xfac * dstress_plus + dstress_minus) * _elasticity_tensor[_qp];
  }
}
