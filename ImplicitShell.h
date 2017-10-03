#include <iostream>
#include <cmath>
#include <Eigen/Geometry>
#include <igl/per_face_normals.h>
#include <igl/per_edge_normals.h>
#include "ArbitraryPrecision.h"

using namespace Eigen;
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Vector6d (*FuncRHS) (int);

void CalcIbar(std::vector<Matrix2m> *abars, std::vector<Matrix2m> *bbars);
MatrixXd CalcMD(MatrixXd GlobalMD);
void Sim();
VectorXd VMass();
DiagonalMatrix<double,Dynamic> InverseMass();
MatrixXd VectoMatrix(VectorXd);
Vector3d edgeV(Vector3d, Vector3d);
Matrix3d rMatrix(Vector3d, Vector3d);
double TotEnergy();
double Energy(double, double, double, int, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d);
Matrix2d QMatrix(double, double, double);
MatrixXd Force();
MatrixXd DelI(Vector3d, Vector3d, Vector3d);
VectorXi NeighborF();
MatrixXi VofEdgeN(); 
Matrix3d CrossM(Vector3d);
Matrix3d dNormVecM(Vector3d);
MatrixXd dEnormM();
Vector2d FoldCircle(double,  double, double);
MatrixXd Enormal(MatrixXd);
Matrix2d I2(Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d);
Matrix2d I1(Vector3d, Vector3d, Vector3d);
void diffusion_prism();
Matrix3d M_Jacobian(int i);
Vector6d UniTopF_prism(int);
VectorXd F_Total();
Matrix6d Mass_prism(int i);
SparseMatrix<double> Mass_Total();
Matrix6d Stiffness_prism(int i);
SparseMatrix<double> Stiffness_Total();
Vector6d PlaneSourceBot(int);

