#ifndef SIMULATIONMESH_H
#define SIMULATIONMESH_H

#include "MaterialParameters.h"
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include "/u/hsiaoyu/Documents/Original_Evouga/ShellOptimization/alglib/src/optimization.h"
#include "ArbitraryPrecision.h"

class RenderingMesh;

class SimulationMesh
{
public:
    SimulationMesh(const Eigen::MatrixX3d &Vbar, const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, double thickness, double youngs, double poisson);

    friend SimulationMesh *parseCirclePacking(const std::string &vertexFile, const std::string &faceFile, const std::string &conformalFactorFile, double scale, double thickness);
    
    const MatrixX3m &getOriginalV() { return originalV; }
    const Eigen::MatrixX3i &getOriginalF() { return F; }

    void faceEnergyDensities(const MatrixX3m &V, VectorXm &densities);

    //void testElasticEnergy();
    //void testElasticEnergy(const std::vector<Matrix2m> &abars, const std::vector<Matrix2m> &bbars);
    void testElasticEnergy(VectorXm *dE, Eigen::SparseMatrix<double> *hEnergy1, Eigen::SparseMatrix<double> *hEnergy2, std::vector<Matrix2m> *acurrent, std::vector<Matrix2m> *bcurrent, const std::vector<Matrix2m> &abars, const std::vector<Matrix2m> &bbars);
    void lineSearch(const MatrixX3m &V, const std::vector<Matrix2m> &targetas, const VectorXm &deltaV, scalar &alpha, scalar &energy);
    void minimizeElasticEnergyNewton(MatrixX3m &V, int substeps, RenderingMesh &rm);
    std::vector<Matrix2m> baras;
    void SetParam(double youngs, double poisson) {params.YoungsModulus = youngs, params.PoissonsRatio = poisson;}

private:
    void buildFaceWings();

    MatrixX3m originalV;
    Eigen::MatrixX3i F;
    Eigen::MatrixX3i faceWings;

    VectorXm faceThicknesses;
    MaterialParameters params;
};

SimulationMesh *parseCirclePacking(const std::string &vertexFile, const std::string &faceFile, const std::string &conformalFactorFile, double scale, double thickness);

#endif
