#include "SimulationMesh.h"
#include <string>
#include <vector>
#include <fstream>
#include <random>
#include "Eigen/Core"
#include "Geometry.h"
#include "ShellEnergy.h"
#include <iostream>
//#include "Eigen/CholmodSupport"
#include "RenderingMesh.h"
#include "ArbitraryPrecision.h"
#include "Eigen/Dense"
#include "CTCD.h"

using namespace std;
using namespace Eigen;
using namespace mpfr;
/*
SimulationMesh *parseCirclePacking(const std::string &vertexFile, const std::string &faceFile, const std::string &conformalFactorFile, double scale, double thickness)
{
    vector<double> xs;
    vector<double> ys;
    ifstream ifs(vertexFile.c_str());
    if (!ifs)
    {
        cout << "Bad vertex file" << endl;
        return NULL;
    }

    while (true)
    {
        double x, y;
        ifs >> x >> y;
        if (!ifs)
            break;
        xs.push_back(x);
        ys.push_back(y);
    }
    int nverts = (int)xs.size();
    MatrixX3d V(nverts, 3);
    V.setRandom();
    V *= 1e-4;
    for (int i = 0; i < nverts; i++)
    {
        V(i, 0) += xs[i];
        V(i, 1) += ys[i];
        V(i, 2) += 0;
    }
    V *= scale;
    vector<int> fs;
    ifstream iff(faceFile.c_str());
    if (!iff)
        return NULL;
    while (true)
    {
        char c;
        int f1, f2, f3;
        iff >> c >> f1 >> f2 >> f3;
        if (!iff)
            break;
        if (c != 'f')
            return NULL;
        fs.push_back(f1 - 1);
        fs.push_back(f2 - 1);
        fs.push_back(f3 - 1);
    }
    int nfaces = (int)fs.size() / 3;
    MatrixX3i F(nfaces, 3);
    for (int i = 0; i < nfaces; i++)
    {
        F(i, 0) = fs[3 * i];
        F(i, 2) = fs[3 * i + 1];
        F(i, 1) = fs[3 * i + 2];
    }

    vector<double> factors;
    ifstream ifc(conformalFactorFile.c_str());
    if (!ifc)
        return NULL;
    for (int i = 0; i < nverts; i++)
    {
        double fac;
        ifc >> fac;
        factors.push_back(fac);
    }
    if (!ifc)
        return NULL;
  */  
    /*V.resize(6, 3);
    V << 0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
        2, 0, 0,
        2, 1, 0;
    F.resize(4, 3);
    F << 0, 1, 2,
        0, 2, 3,
        1, 4, 2,
        2, 4, 5;
    nfaces = 4;*/
/*
    SimulationMesh *result = new SimulationMesh(V, F, thickness);

    for (int i = 0; i < nfaces; i++)
    {
        double avfac = 0;
        for (int j = 0; j < 3; j++)
        {
            avfac += factors[F(i, j)];
        }
        avfac /= 3.0;
        result->baras[i] *= avfac;
    }

    return result;
}
*/
SimulationMesh::SimulationMesh(const Eigen::MatrixX3d & Vbar, const Eigen::MatrixX3d & V, const Eigen::MatrixX3i & Fin, double thickness)
{
    int nverts = (int)V.rows();
    int nfaces = (int)Fin.rows();
    originalV.resize(nverts, 3);
    for (int i = 0; i < nverts; i++)
        for (int j = 0; j < 3; j++)
            originalV(i, j) = V(i, j);
    F.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++)
        for (int j = 0; j < 3; j++)
            F(i, j) = Fin(i, j);
    baras.resize(nfaces);
    faceThicknesses.resize(nfaces);
    for (int i = 0; i < nfaces; i++)
    {
        Vector3i face = F.row(i);
        baras[i] = firstFundamentalForm(Vbar.row(face[0]), Vbar.row(face[1]), Vbar.row(face[2]));
        faceThicknesses[i] = thickness;
    }

    buildFaceWings();
}

/*  // The original copy of the code
SimulationMesh::SimulationMesh(const Eigen::MatrixX3d & V, const Eigen::MatrixX3i & Fin, double thickness)
{
    int nverts = (int)V.rows();
    int nfaces = (int)Fin.rows();
    originalV.resize(nverts, 3);
    for (int i = 0; i < nverts; i++)
        for (int j = 0; j < 3; j++)
            originalV(i, j) = V(i, j);
    F.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++)
        for (int j = 0; j < 3; j++)
            F(i, j) = Fin(i, j);
    baras.resize(nfaces);
    faceThicknesses.resize(nfaces);
    for (int i = 0; i < nfaces; i++)
    {
        Vector3i face = F.row(i);
        baras[i] = firstFundamentalForm(originalV.row(face[0]), originalV.row(face[1]), originalV.row(face[2]));
        faceThicknesses[i] = thickness;
    }

    buildFaceWings();
}
*/

void SimulationMesh::buildFaceWings()
{
    int nfaces = (int)F.rows();
    faceWings.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int result = -1;
            int p1 = F(i, (j + 1) % 3);
            int p2 = F(i, (j + 2) % 3);
            for (int k = 0; k < nfaces; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    if (F(k, (l + 1) % 3) == p2 && F(k, (l + 2) % 3) == p1)
                    {
                        result = F(k, l);
                    }
                }
            }
            faceWings(i, j) = result;
        }
    }
}

//void SimulationMesh::testElasticEnergy()
//void SimulationMesh::testElasticEnergy(const std::vector<Matrix2m> &abars, const std::vector<Matrix2m> &bbars)
void SimulationMesh::testElasticEnergy(VectorXm *dE, Eigen::SparseMatrix<double> *hEnergy1, Eigen::SparseMatrix<double> *hEnergy2, const std::vector<Matrix2m> &abars, const std::vector<Matrix2m> &bbars)
{
    int nDOF = originalV.rows() *3;
    dE->resize(nDOF);
    dE->setZero();
    hEnergy1->resize(nDOF, nDOF );
    hEnergy1->setZero();
    hEnergy2->resize(nDOF, nDOF );
    hEnergy2->setZero();
    std::vector<Eigen::Triplet<scalar> > hEnergyList1, hEnergyList2;
    scalar result = shellEnergy(originalV, F, faceWings, baras, bbars, faceThicknesses, params, dE, &hEnergyList1, &hEnergyList2, NULL);
    hEnergy1->setFromTriplets(hEnergyList1.begin(), hEnergyList1.end());
    hEnergy2->setFromTriplets(hEnergyList2.begin(), hEnergyList2.end());
    //cout << "Force: " << endl << *dE << endl;
/*
    int nverts = (int)originalV.rows();
    for (int i = 0; i < 3*nverts; i++)
    {
        MatrixX3m newV = originalV;
        newV(i/3, i%3) += 1e-8;
        scalar newe = shellEnergy(newV, F, faceWings, baras, faceThicknesses, params, NULL, NULL, NULL, NULL);
        cout << (newe - result) / 1e-8 << " vs " << dE[i] << endl;        
    }
*/
}

//void SimulationMesh::lineSearch(const MatrixX3m &V, const std::vector<Matrix2m> &targetas, const VectorXm &deltaV, scalar &t, scalar &energy)
//{
//    t = 1.0;
//    scalar c1 = 0.8;
//    scalar c2 = 0.9;
//    scalar alpha = 0;
//    scalar infinity = 1e6;
//    scalar beta = infinity;
//
//    VectorXm dE;
//    VectorXm newdE;
//    dE.setZero();
//    scalar orig = shellEnergy(V, F, faceWings, targetas, faceThicknesses, params, &dE, NULL, NULL,NULL);
//    scalar deriv = dE.dot(deltaV);
//    assert(deriv < 0);
//    int nverts = (int)V.rows();
//    while (true)
//    {
//        MatrixX3m testV = V;
//        for (int i = 0; i < nverts; i++)
//        {
//            testV.row(i) += t*deltaV.segment<3>(3 * i);
//        }
//        newdE.setZero();
//        scalar newenergy = shellEnergy(testV, F, faceWings, targetas, faceThicknesses, params, &newdE, NULL, NULL, NULL);
//        if (isnan(newenergy) || newenergy > orig + t*deriv*c1)
//        {
//            beta = t;
//            t = 0.5*(alpha + beta);
//        }
//        else if (newdE.dot(deltaV) < c2*deriv)
//        {
//            alpha = t;
//            if (beta == infinity)
//            {
//                t = 2 * alpha;
//            }
//            else
//            {
//                t = 0.5*(alpha + beta);
//            }
//
//            if (beta - alpha < 1e-8)
//            {
//                energy = newenergy;
//                return;
//            }
//        }
//        else
//        {
//            energy = newenergy;
//            return;
//        }
//    }
//}
//
//double powerIteration(const Eigen::SparseMatrix<double> &M, double tol, double bias, Eigen::VectorXd &v)
//{
//    int ndofs = (int)M.rows();
//    VectorXd seed(ndofs);
//    seed.setRandom();
//    seed = seed / seed.norm();
//    double lambda = 0;
//    while (true)
//    {
//        double newlambda = seed.norm();
//        if (fabs(lambda - newlambda) < tol)
//        {
//            v = seed;
//            return seed.dot(M*seed + bias*seed) / newlambda / newlambda;
//        }
//        //cout << fabs(lambda - newlambda) << endl;
//        lambda = newlambda;
//        seed /= lambda;
//        seed = M*seed + bias*seed;
//    }
//}
//
//void SimulationMesh::minimizeElasticEnergyNewton(MatrixX3m &V, int substeps, RenderingMesh &rm)
//{
//    int nverts = (int)V.rows();
//    int nfaces = (int)F.rows();
//    for (int i = 0; i < substeps; i++)
//    {
//        cout << "Starting substep " << i << endl;       
//        vector<Matrix2m> targetas(nfaces);
//        for (int j = 0; j < nfaces; j++)
//        {
//            Vector3i face = F.row(j);
//            Matrix2m origa = firstFundamentalForm(V.row(face[0]), V.row(face[1]), V.row(face[2]));
//            Matrix2m term1 = baras[j] * scalar(i + 1) * scalar(i + 1);
//            Matrix2m term2 = origa * scalar(substeps*substeps - (i + 1)*(i + 1));
//            targetas[j] = (term1 + term2) / scalar(substeps*substeps);
//        }
//        
//        scalar dEnorm = numeric_limits<double>::infinity();
//        int iter = 0;
//        SimplicialLDLT < SparseMatrix < scalar > > solver;
//        
//        scalar reg = 1.0;
//        while (true)
//        {
//            iter++;
//            VectorXm dE(3*nverts);
//            dE.setZero();
//            vector<Eigen::Triplet<scalar> > hEexact;
//            vector<Eigen::Triplet<scalar> > hEinexact;
//            VectorXm triEnergies(nfaces);
//            triEnergies.setZero();
//            cout << "Computing Hessian" << endl;            
//            scalar result = shellEnergy(V, F, faceWings, targetas, faceThicknesses, params, &dE, &hEexact, &hEinexact, &triEnergies);            
//
//            dEnorm = dE.norm();
//            cout << "Iter " << iter << " energy: " << result << " grad " << dEnorm << endl;
//
//            SparseMatrix <scalar > hessexact(3 * nverts, 3 * nverts);            
//            hessexact.setFromTriplets(hEexact.begin(), hEexact.end());   
//            
//            SparseMatrix < scalar > hessinexact(3 * nverts, 3 * nverts);
//            hessinexact.setFromTriplets(hEinexact.begin(), hEinexact.end());
//            
//            VectorXm rhs = -dE;
//            
//            cout << "Factorizing Hessian, regularization = " << reg << endl;
//            solver.setShift(reg);
//            solver.compute(hessexact);
//            
//            while (true)
//            {
//                if (solver.info() == Success)
//                {
//                    bool hasneg = false;
//                    for (int i = 0; i < 3 * nverts; i++)
//                        if (solver.vectorD()[i] < -1e-4)
//                            hasneg = true;
//                    if (!hasneg)
//                        break;
//                }
//                
//                reg *= 2.0;
//                cout << "Factorizing Hessian, regularization = " << reg << endl;
//                solver.setShift(reg);
//                solver.compute(hessexact);
//            }
//
//            if(reg > 1e-8)
//                reg /= 2.0;
//
//            bool exact = false;
//            if (solver.info() == Success)
//            {                
//                cout << "Solving, exact Hessian" << endl;
//                bool hasneg = false;
//                for (int i = 0; i < 3 * nverts; i++)
//                    if (solver.vectorD()[i] < -1e-4)
//                        hasneg = true;
//                if(!hasneg)
//                    exact = true;
//            }
//            if(!exact)
//            {
//                cout << "Solving, inexact Hessian" << endl;
//                solver.compute(hessinexact);
//                for (int i = 0; i<3 * nverts; i++)
//                    if (solver.vectorD()[i] < 0)
//                        cout << solver.vectorD()[i] << endl;
//            }
//
//            if (solver.info() != Success)
//            {
//                cout << "Cannot factorize inexact Hessian, error " << solver.info() << endl;
//                MatrixXd dense(hessinexact);
//                cout << "build dense matrix" << endl;
//                Eigen::SelfAdjointEigenSolver<MatrixXd> sol(dense);
//                cout << sol.eigenvalues().transpose() << endl;
//                while (true);
//            }
//
//       
//             
//            VectorXm newtondeltaV = solver.solve(rhs);
//            
//            //cout << newtondeltaV.transpose() << endl;
//
//            std::cout << "Descent dir magnitude " << newtondeltaV.norm() << endl;
//            std::cout << "Performing line search" << endl;
//
//            scalar newtonalpha = 0;
//            scalar newtonenergy = 0;
//            if (newtondeltaV.dot(dE) > 0)
//                newtonenergy = numeric_limits<double>::infinity();
//            else
//                lineSearch(V, targetas, newtondeltaV, newtonalpha, newtonenergy);
//            
//            scalar descentalpha = 0;
//            scalar descentenergy = 0;
//            
//            //lineSearch(V, targetas, -dE, descentalpha, descentenergy);
//
//            cout << "Newton step " << newtonalpha << " energy diff " << result - newtonenergy << endl;
//            //cout << "Descent step " << descentalpha << " energy diff " << result - descentenergy << endl;
//
//            scalar alpha;
//            VectorXm deltaV;
//            /*if (descentenergy < newtonenergy)
//            {
//                alpha = descentalpha;
//                deltaV = -dE;
//            }
//            else*/
//            {
//                alpha = newtonalpha;
//                deltaV = newtondeltaV;
//            }            
//
//            rm.visMutex.lock();
//            rm.V.resize(V.rows(), 3);
//            for (int i = 0; i < (int)V.rows(); i++)
//                for (int j = 0; j < 3; j++)
//                    rm.V(i, j) = double(V(i, j));
//            rm.faceVals.resize(triEnergies.size());
//            for (int i = 0; i<(int)triEnergies.size(); i++)
//                rm.faceVals[i] = double(triEnergies[i]);
//            for (int i = 0; i < nfaces; i++)
//            {
//                double det = baras[i].determinant();
//                double area = 0.5*std::sqrt(det);
//                rm.faceVals[i] /= area;
//            }
//            rm.gradientDir.resize(3 * nverts);
//            for (int i = 0; i < 3 * nverts; i++)
//                rm.gradientDir[i] = double(dE[i]);
//            rm.descentDir.resize(3 * nverts);
//            for (int i = 0; i < 3 * nverts; i++)
//                rm.descentDir[i] = double(deltaV[i]);
//            rm.visMutex.unlock();
//            
//            for (int i = 0; i < nverts; i++)
//            {
//                V.row(i) += alpha*deltaV.segment<3>(3 * i);
//            }
//
//            if(exact && reg <= 1e-8 && newtondeltaV.norm()/result < 1e-6)
//                break;
//        }
//        rm.writeMesh(i);
//    }
//}
//
//void SimulationMesh::faceEnergyDensities(const MatrixX3m &V, VectorXm &densities)
//{
//    int nfaces = F.rows();
//    densities.resize(nfaces);
//    densities.setZero();
//    shellEnergy(V, F, faceWings, baras, faceThicknesses, params, NULL, NULL, NULL, &densities);
//    for (int i = 0; i < nfaces; i++)
//    {
//        double area = 0.5*std::sqrt(baras[i].determinant());
//        densities[i] /= area;
//    }
//}
