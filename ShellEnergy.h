#ifndef SHELLENERGY_H
#define SHELLENERGY_H

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "ArbitraryPrecision.h"
#include <vector>

struct MaterialParameters;

scalar triangleStretchingEnergy(
    const Vector3m &p1,
    const Vector3m &p2,
    const Vector3m &p3,
    const Matrix2m &abar,
    scalar h,
    const MaterialParameters &params,
    Eigen::Matrix<scalar, 9, 1> *dEnergy,
    Eigen::Matrix<scalar, 9, 9> *hEnergyExact,
    Eigen::Matrix<scalar, 9, 9> *hEnergyInexact);

scalar triangleBendingEnergy(
    const Vector3m &p0,
    const Vector3m &p1,
    const Vector3m &p2,
    const Vector3m *q0,
    const Vector3m *q1,
    const Vector3m *q2,
    const Matrix2m &abar,
    const Matrix2m &bbar,
    scalar h,
    const MaterialParameters &params,
    Matrix2m *a,
    Matrix2m *b,
    Eigen::Matrix<scalar, 9, 1> *dEnergy,
    Eigen::Matrix<scalar, 9, 9> *hEnergyInexact);


scalar shellEnergy(
    const MatrixX3m &V,
    const Eigen::MatrixX3i &F,
    const Eigen::MatrixX3i &faceWings,
    const std::vector<Matrix2m> &abars,
    const std::vector<Matrix2m> &bbars,
    const VectorXm &faceThicknesses,
    const MaterialParameters &params,
    std::vector<Matrix2m> *acurrent,
    std::vector<Matrix2m> *bcurrent,
    VectorXm *dEnergy,
    std::vector<Eigen::Triplet<scalar> > *hEnergyExact,
    std::vector<Eigen::Triplet<scalar> > *hEnergyInexact,
    VectorXm *triangleEnergies);

scalar poseEnergy(
    const MatrixX3m &V,
    const Eigen::MatrixX3i &F,
    const MatrixX3m &origV,
    VectorXm *dEnergy,
    std::vector<Eigen::Triplet<scalar> > *hEnergyExact,
    std::vector<Eigen::Triplet<scalar> > *hEnergyInexact);

#endif
