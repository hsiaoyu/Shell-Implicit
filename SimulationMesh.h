#ifndef SIMULATIONMESH_H
#define SIMULATIONMESH_H

#include <Eigen/Sparse>
#include "MaterialParameters.h"
#include <vector>
#include "/u/hsiaoyu/Documents/evouga/ShellOptimization/alglib/src/optimization.h"
#include "ArbitraryPrecision.h"

class RenderingMesh;

class SimulationMesh
{
public:
    SimulationMesh(Eigen::MatrixX3d Vbar, Eigen::MatrixX3d V, Eigen::MatrixX3i F, double thickness);
    friend SimulationMesh *parseCirclePacking(const std::string &vertexFile, const std::string &faceFile, const std::string &conformalFactorFile, double scale, double thickness);
    
    const MatrixX3m &getOriginalV() { return originalV; }
    const Eigen::MatrixX3i &getOriginalF() { return F; }

    void faceEnergyDensities(const MatrixX3m &V, VectorXm &densities);

    void testElasticEnergy(VectorXm *dE, Eigen::SparseMatrix<double> *hEnergy1, Eigen::SparseMatrix<double> *hEnergy2);
    void lineSearch(const MatrixX3m &V, const std::vector<Matrix2m> &targetas, const VectorXm &deltaV, scalar &alpha, scalar &energy);
    void minimizeElasticEnergyNewton(MatrixX3m &V, int substeps, RenderingMesh &rm);
    std::vector<Matrix2m> baras;

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
