// The file generate the animation of the shell with moisture diffusing throughout the material
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_edge_normals.h>
#include "ImplicitShell.h"
#include "SimulationMesh.h"
#include "RenderingMesh.h"
#include <thread>
#include <igl/writeOBJ.h>
#include <vector>

using namespace std;
using namespace Eigen;

Matrix2d CanonM1, CanonM2, CanonM3;
double area=0.5, constE, nu, t, delt, ShrinkCoeffMD, ShrinkCoeffCD, rho, DiffusionCoeff, Dissipation, cdamp1, cdamp2, sourcetime, SourceZ, LameFirst, LameSecond, ImpTolerance; // area of the canonical triangle
MatrixXd velocity, V, Vbar, ENbar, Ibartot, MachineDirection, IAold, IBold;// ENbar is the edge normal of Vbar
MatrixXi F, EdgeNV; // Eash row of EdgeNV(3*NTri,4) stores the index off the four vertices that are related to one edge
VectorXi EdgeF;// EdgeF(3*NTri) stores the index of the  adjecent face of the edge other than the face that is indicated by the vector index 
int NNode, gflag, tcount, storefreq, InitialNum, SourceType; // NNode is the total # of nodes
Vector3d V_ini1, V_ini2; // Vini1, Vini2 are the fixed points of V when enabled gravity
VectorXd Moisture_v, RHS_F, Mass, Pos, Vel, Moist_source;
SparseMatrix<double> Mtot_Mass, Mtot_Stiffness;
MatrixXd Vtime; // Vtime stores the vertices at each time 
string objname, moistname;
bool DampingForce_Enabled, Implicit, forceMD=0;
FuncRHS funcRHS[1]={UniTopF_prism};
DiagonalMatrix<double,Dynamic> MassInv, MassVec;
vector<Matrix2m> a_prev, b_prev;

ofstream myfile;
int main(int argc, char* argv[]){
        CanonM3 << 1, 0, 0, 0;
        CanonM2 << 1, 1, 1 ,1;
        CanonM1 << 0, 0, 0, 1;
        int i, j, tmp, NTri;
        Vector3d Ei, Ef;
    //    ifstream Infile("burning_diffusion.txt");
        ifstream Infile(argv[1]);
        string str;
        getline(Infile,str);
        const char *cstr = str.c_str();
        MatrixXd GlobalMD;
        igl::readDMAT(cstr, GlobalMD);
        //---------------------------------------------
        getline(Infile,str); 
        const char *cstr1 = str.c_str();
        igl::readOBJ(cstr1,Vbar,F);
        getline(Infile,str); 
        const char *cstr2 = str.c_str();
	igl::readOBJ(cstr2,V,F);
        getline(Infile,str); 
        const char *cstr3 = str.c_str();
        // Moisture_v(i) stores the moist of vertice i and Moisture_v(i+NNode) stores the value for the corresponding lower surface
        igl::readDMAT(cstr3, Moisture_v);
        Moist_source = Moisture_v; // If a source exist the moisture to be added at each time step would be the same as the initial source
        getline(Infile,objname); 
        getline(Infile,moistname); 
        Infile >> constE;
        Infile >> nu;
        Infile >> t;
        Infile >> rho;
        Infile >> delt;
	// include gravity with gflag = 1
        Infile >> Implicit;
        Infile >> ImpTolerance;
        Infile >> gflag;
        Infile >> ShrinkCoeffMD;
        Infile >> ShrinkCoeffCD;
	Infile >> DiffusionCoeff;
	Infile >> DampingForce_Enabled;
	Infile >> cdamp1;
	Infile >> cdamp2;
	Infile >> storefreq; // The frequency of storing V to Vtime for animation
	Infile >> InitialNum; // starting file number for recording the animation 
	// 0-No source, 1-Const Dissipation, 2- Source at certain Z-plane
	Infile >> SourceType; // What kind of source for moisture 
	Infile >> Dissipation;// The (negative)amount of moisture added to the material per second 
	Infile >> sourcetime; // Total amount of time the source would act on the shell, for time longer than sourcetime source would be set to 0 
	Infile >> SourceZ; // The Z-coordinate of the source 
	Infile >> forceMD; // TODO: remove this 

	LameFirst = constE*nu/(1+nu)/(1-2*nu);
	LameSecond = constE/2/(1+nu);
	
	NTri=F.rows();
	NNode=V.rows();	
        Pos.resize(3*NNode);
        Vel.setZero(3*NNode);

	// Initialize position and velocity in vector form for implicit calculation
        for (int i=0; i<NNode; i++){
	    for (int j=0; j<3; j++){
		Pos(3*i+j)=V(i,j);
	    }
	}
        Mass = VMass();
	MassInv = InverseMass();
	MassVec = MasstoMatrix();
	//TODO: Implement the correct MD calc for 3D
	if (forceMD){
	  MachineDirection.resize(2*NTri,2);
	  for(int i=0; i<NTri; i++){
	     MachineDirection(i,0) = 1;
	     MachineDirection(i,1) = 0;
	     MachineDirection(i+NTri,0) = 0;
	     MachineDirection(i+NTri,1) = 1;
          } 
	}
	else{
          MachineDirection = CalcMD(GlobalMD);
	}
        EdgeF=NeighborF(); 
        EdgeNV=VofEdgeN();
        MatrixXd FNbar(NTri,3);
	igl::per_face_normals(Vbar, F, FNbar);
        ENbar=Enormal(FNbar);
        velocity.setZero(NNode,3);
        RHS_F.resize(2*NNode,1);
        RHS_F = F_Total();
	Mtot_Mass=Mass_Total();
    	Mtot_Stiffness=Stiffness_Total();
    	Mtot_Stiffness*=DiffusionCoeff;
        IAold.resize(2*NTri,2);
        IBold.resize(2*NTri,2);
	tcount =0; // Total time step
	while(true){
	  if ((tcount%storefreq)==0){
	  	stringstream FileName;
	  	FileName << objname << tcount/storefreq+InitialNum << ".obj"; 
	  	igl::writeOBJ(FileName.str(),V,F);
	  	FileName.str("");
	  	FileName << moistname << tcount/storefreq+InitialNum << ".dmat"; 
	  	igl::writeDMAT(FileName.str(),Moisture_v,1);
	  }
	  Sim();
  	  tcount++;
	}

	return 0;
}
      
//CalcIAbar calculates IAbar and IBbar based on the moisture distribution assuming that only upper surface changes and the lower surface remains intact
void CalcIbar(vector<Matrix2m> *abars, vector<Matrix2m> *bbars){
	int NTri = F.rows();
	MatrixXd IAtot(2*NTri,2), IBtot(2*NTri,2), Itot(4*NTri,2),FN(NTri,3), EN(3*NTri,3);
	Matrix2m IAbarnew, IBbarnew;
	//abars.resize(NTri);
	//bbars.resize(NTri);
//	igl::per_face_normals(V, F, FN);
//	EN=Enormal(FN);
	for (int i=0; i<NTri; i++){
		Matrix2d IA, IB, tmp1, tmp2, U, Utransform, MShrink;
		double c1, c2, c3, ctemp, MoistureLevel;
		Vector2d MD, CD;
		// Moisture level on a face is defined by averaging the moisture on the three vertices
		MoistureLevel=(double)(Moisture_v(F(i,0))+Moisture_v(F(i,1))+Moisture_v(F(i,2)))/3;
		IB=I2(Vbar.row(F(i,0)),Vbar.row(F(i,1)),Vbar.row(F(i,2)),ENbar.row(3*i),ENbar.row(3*i+1),ENbar.row(3*i+2));
		IA=I1(Vbar.row(F(i,0)),Vbar.row(F(i,1)),Vbar.row(F(i,2)));
		//IB=I2(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),EN.row(3*i),EN.row(3*i+1),EN.row(3*i+2));
		//IA=I1(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
//	Remove MD CD for now
 		tmp1=IA+t*IB;
		tmp2=IA-t*IB;
		//Canonical e1=(1,0), e2=(-1,1), MachineDirection is the coefficient of the globe MD on the barycentric coord of this triangle.
		MD<<MachineDirection(i,0)-MachineDirection(i,1),MachineDirection(i,1);
		CD<<MachineDirection(i+NTri,0)-MachineDirection(i+NTri,1),MachineDirection(i+NTri,1);
		//----------------------Method 1----------------------------------------
  		ctemp=MD.transpose()*tmp1*MD;
		c1=(1.0-ShrinkCoeffMD*MoistureLevel)*(1.0-ShrinkCoeffMD*MoistureLevel)*ctemp;
		ctemp=CD.transpose()*tmp1*MD;
		c2=(1.0-ShrinkCoeffCD*MoistureLevel)*(1.0-ShrinkCoeffMD*MoistureLevel)*ctemp;
		ctemp=CD.transpose()*tmp1*CD;
		c3=(1.0-ShrinkCoeffCD*MoistureLevel)*(1.0-ShrinkCoeffCD*MoistureLevel)*ctemp;
		U << MD(0), CD(0), MD(1), CD(1);
		Utransform=U.inverse(); 
		//Utransform *= 1/(MD(0)*CD(1)-MD(1)*CD(0));
		MShrink << c1, c2, c2, c3;
		tmp1= Utransform.transpose()*MShrink*Utransform;
		// Calculate I for the lower surface
                ctemp=MD.transpose()*tmp2*MD;
		MoistureLevel=(double)(Moisture_v(F(i,0)+NNode)+Moisture_v(F(i,1)+NNode)+Moisture_v(F(i,2)+NNode))/3;
		c1=(1.0-ShrinkCoeffMD*MoistureLevel)*(1.0-ShrinkCoeffMD*MoistureLevel)*ctemp;
		ctemp=CD.transpose()*tmp2*MD;
		c2=(1.0-ShrinkCoeffCD*MoistureLevel)*(1.0-ShrinkCoeffMD*MoistureLevel)*ctemp;
		ctemp=CD.transpose()*tmp2*CD;
		c3=(1.0-ShrinkCoeffCD*MoistureLevel)*(1.0-ShrinkCoeffCD*MoistureLevel)*ctemp;
		MShrink << c1, c2, c2, c3;
		tmp2= Utransform.transpose()*MShrink*Utransform;
		IAbarnew = (tmp1+tmp2)/2;
		IBbarnew = (tmp1-tmp2)/(2*t);
		//--------------------------------------------------------------------------
		
		/*//-----------Method 2--------------------------------------
		Matrix2d T;
 		tmp1=IA+t*IB;
		tmp2=IA-t*IB;
		// Calculate top layer
		MoistureLevel=(double)(Moisture_v(F(i,0))+Moisture_v(F(i,1))+Moisture_v(F(i,2)))/3;
		c1=(1.0-ShrinkCoeffMD*MoistureLevel);
		c2=(1.0-ShrinkCoeffCD*MoistureLevel);
		U << MD(0), CD(0), MD(1), CD(1);
		MShrink << c1, 0, 0, c2;
		Utransform=U.inverse();
		T = U*MShrink*Utransform; 
		tmp1 = T.transpose()*tmp1*T;
		// Calculate bottom layer
		MoistureLevel=(double)(Moisture_v(F(i,0)+NNode)+Moisture_v(F(i,1)+NNode)+Moisture_v(F(i,2)+NNode))/3;
		c1=(1.0-ShrinkCoeffMD*MoistureLevel);
		c2=(1.0-ShrinkCoeffCD*MoistureLevel);
		U << MD(0), CD(0), MD(1), CD(1);
		MShrink << c1, 0, 0, c2;
		Utransform=U.inverse();
		T = U*MShrink*Utransform; 
		tmp2 = T.transpose()*tmp2*T;
		IAbarnew = (tmp1+tmp2)/2;
		IBbarnew = (tmp1-tmp2)/(2*t);
		//-----------------------------------------------------------
		*/
		IAtot.block(2*i,0,2,2)=IAbarnew;
		IBtot.block(2*i,0,2,2)=IBbarnew;
	//	abars[i]=IAbarnew;	
	//	bbars[i]=IBbarnew;
		abars->push_back(IAbarnew);	
		bbars->push_back(IBbarnew);	
		
		/* // Sometimes with a bad mesh there would be a slight discrepancy between IAbarnew and IA even with moisture 
 * 		// difference is zero. Below is the test code for the discrepancy.  
		if ((IAbarnew-IA).lpNorm<Infinity>() > 1e-16){
		cout << i << endl;
		cout << F(i,0) << " " << F(i,1) << " " << F(i,2) << endl;
		cout <<"IA" << endl << IAbarnew -IA << endl;
		}
		if ((IBbarnew-IB).lpNorm<Infinity>() > 1e-16){
		cout << i << endl;
		cout << F(i,0) << " " << F(i,1) << " " << F(i,2) << endl;
		cout <<"IB" << endl << IBbarnew -IB << endl;
		}
		*/
		//IAtot.block(2*i,0,2,2)=IA;
		//IBtot.block(2*i,0,2,2)=IB;
		
	}

	Itot << IAtot, IBtot;
	Ibartot = Itot;
}

MatrixXd CalcMD3d(Vector3d GMD){
   int NTri = F.rows();
   MatrixXd BaryMD(NTri,2);
   Vector3d cb1, cb2, cb3;
   cb1 = V.row(F(0,1))-V.row(F(0,0));
   cb2 = V.row(F(0,2))-V.row(F(0,1));
   cb3 = cb1.cross(cb2);
   Matrix3d B;
   B.col(0) = cb1;  
   B.col(1) = cb2;  
   B.col(2) = cb3;  
   Vector3d tmp = B.inverse() * GMD;
   BaryMD(0,0) = tmp(0); 
   BaryMD(0,1) = tmp(1);
   int CalcedFace[NTri];
   for (int i=0; i<NTri; i++){
       CalcedFace[i]=0;
   }
   CalcNextMD(&BaryMD,0,CalcedFace);  
}

void CalcNextMD(MatrixXd *BaryMD, int FaceID, int CalcedFace[]){
  Vector3d curMD;
  curMD << BaryMD->row(FaceID), 0;
  for (int i=0; i<3 ; i++){
      int curF = EdgeF(3*FaceID+i); 
      if( curF > 0 && CalcedFace[curF]==0){
        Vector3d tmp;
        tmp = CalcNeighborMD(curMD, FaceID, i);
	BaryMD->row(curF) << tmp(0), tmp(1);
        CalcedFace[curF]=1;
        CalcNextMD(BaryMD, curF, CalcedFace);
      }
   }
}


Vector3d CalcNeighborMD(Vector3d curMD, int currentF, int EdgeID){
   Vector3d cb1, cb2, cb3, edge, ce2;
   // cb1, cb2 are the barycentric vector of the current face
   cb1 = V.row(F(currentF,1))-V.row(F(currentF,0));
   cb2 = V.row(F(currentF,2))-V.row(F(currentF,1));
   cb3 = cb1.cross(cb2);
   if(EdgeID==0){
     edge = cb1;
   }
   else if(EdgeID==1){
     edge = cb2;
   }
   else if(EdgeID==2){
     edge = V.row(F(currentF,0))-V.row(F(currentF,2));
   }
   Matrix3d rot;
   rot = AngleAxisd(0.5*M_PI,cb3);
   // ce2 is perpendicular to the shared edge and lay on the current face
   ce2 = rot * edge;
   
   int neighborF = EdgeF(3*currentF+EdgeID);
   Vector3d nb1, nb2, nb3, ne2;
   nb1 = V.row(F(neighborF,1))-V.row(F(neighborF,0));
   nb2 = V.row(F(neighborF,2))-V.row(F(neighborF,1));
   nb3 = nb1.cross(nb2);
   rot = AngleAxisd(0.5*M_PI,nb3);
   // ne2 is perpendicular to the shared edge and lay on the neighboring face
   ne2 = rot * edge;
   
   Matrix3d Bcur, Ecur, Enei, Bnei;
   Bcur.col(0) = cb1;  
   Bcur.col(1) = cb2;  
   Bcur.col(2) = cb3;  
   Ecur.col(0) = edge;  
   Ecur.col(1) = ce2;  
   Ecur.col(2) = cb3;  
   Bnei.col(0) = nb1;  
   Bnei.col(1) = nb2;  
   Bnei.col(2) = nb3;  
   Enei.col(0) = edge;  
   Enei.col(1) = ne2;  
   Enei.col(2) = nb3;

   Vector3d nMD;
   nMD = Bnei.inverse() * Enei * Ecur.inverse() * Bcur * curMD; 
   
   return nMD;
}
//CalcMD calculates GlobalMD in barycentric coordinate in each face
MatrixXd CalcMD(MatrixXd GlobalMD){
	Vector3d e1, e2;
	Matrix2d A, tmp;
	MatrixXd MD(2*F.rows(),2);
        Vector2d GMD,GCD;
//	GMD is defined as follow if MD is uniform throughout the paper
	GMD << GlobalMD(0,0), GlobalMD(0,1);
        Matrix2d Rot;
	Rot << 0,-1,1,0;
        GCD = Rot*GMD;
	for (int i=0 ; i<F.rows(); i++){
        /*	// If MD is not uniform throughout the paper use the below definition for GMD
		GMD << GlobalMD(i,0), GlobalMD(i,1);
	*/
		e1 = V.row(F(i,1))-V.row(F(i,0));
        	e2 = V.row(F(i,2))-V.row(F(i,1));
		
		tmp << e1(0), e2(0), e1(1), e2(1);
		A = tmp.inverse();
	//	A << e2(1), -e2(0), -e1(1), e1(0);
	//	A *= 1/(e1(0)*e2(1)-e1(1)*e2(0));	
	//	cout << "Inv" << endl <<  Inv << endl << "A" << endl << A << endl;

		Vector2d tmp1, tmp2;
		tmp1 = A*GMD;
                tmp2 = A*GCD;
		// MD is normalized so that IAbar would be easier to calculate by unitary transform
		MD.row(i) = tmp1.transpose()/tmp1.norm();
		MD.row(i+F.rows()) = tmp2.transpose()/tmp2.norm();
//		cout << MD(i,0)*e1+MD(i,1)*e2 << endl << "GMD" << endl <<  GMD << " " << GlobalMD(0,2) << endl;
	}
	return MD;
}

bool pre_draw(igl::viewer::Viewer& viewer){
	if ((tcount%storefreq)==0){
		stringstream FileName;
		FileName << objname << tcount/storefreq+InitialNum << ".obj"; 
		igl::writeOBJ(FileName.str(),V,F);
		FileName.str("");
		FileName << moistname << tcount/storefreq+InitialNum << ".dmat"; 
		igl::writeDMAT(FileName.str(),Moisture_v,1);
	}

	Sim();
	tcount++;

        viewer.data.clear();
    	viewer.data.set_mesh(V, F);
    	viewer.core.align_camera_center(V,F);
	return false;
}


void Sim(){
        VectorXd gravity(NNode);
        gravity.setConstant(100*delt);
	vector<Matrix2m> abars, bbars;
        CalcIbar(&abars, &bbars);
	
	//----------------------Update with Implicit Euler---------------------------------------
	if(Implicit == 1){
 	   
           // Declare a sparse identity matrix of the right size for convenience
           const SparseMatrix<double> identityMatrix = MatrixXd::Identity(3*NNode,3*NNode).sparseView();


           // start with dummy value for residual that's larger than tolerance
           double residual = 1000 + ImpTolerance;

           // Used to store the result of f(guess)
           VectorXd fguess(3*NNode);

           // A biconjugate gradients solver
           //SparseQR<SparseMatrix<double>, COLAMDOrdering<int> >  solver;
           SimplicialLDLT<SparseMatrix<double> >  solver;
           // Solve with Newton's Method
	   SimulationMesh *sm;
           
	   //Use explicit Euler as the initial guess
           VectorXd guess = VectorXd(Pos);
           guess = guess + delt *  Vel;
	   MatrixXd Vtemp;
	   int count = 0;
           while(residual > ImpTolerance){
             VectorXm dE(3*NNode), dE_Damp(3*NNode);
             SparseMatrix<double> hEnergy1(3*NNode,3*NNode), hEnergy2(3*NNode,3*NNode),hDamp1(3*NNode,3*NNode), hDamp2(3*NNode,3*NNode);
	     vector<Matrix2m> atemp, btemp;
 	     Vtemp = VectoMatrix(guess);
             sm = new SimulationMesh(Vbar, Vtemp, F, t, constE, nu);
             sm->testElasticEnergy(&dE, &hEnergy1, &hEnergy2, &atemp, &btemp, abars, bbars);
	     //Calculate damping force
	     //TODO: add comment for only using cdamp1 in Implicit calculation
             if(DampingForce_Enabled && tcount >0){
                sm->testElasticEnergy(&dE_Damp, &hDamp1, &hDamp2, &atemp, &btemp, a_prev, b_prev);
                dE += cdamp1 / delt * dE_Damp;
                hEnergy1 += cdamp1 / delt * hDamp1;
                hEnergy2 += cdamp1 / delt * hDamp2;
             }
	     //Calculate force differential
             SparseMatrix<double> dForces;
             dForces = MassVec * identityMatrix + delt * delt * hEnergy1;
             //Calculate f(guess)
             fguess = MassVec * guess - MassVec * Pos - MassVec * delt *  Vel + delt * delt * dE;

             // Solve dat!
             solver.compute(dForces);
             VectorXd upd = solver.solve(fguess);

             guess -= upd;

             // Calculate residual
             residual = fguess.norm();
             
             if (tcount%100==0 && count%50==0){
                cout <<" tcount = " << tcount << "count = " << count << endl << residual << endl;
	     }
	     count++;
           
             delete sm;
	   }
           Pos = guess;
	   V = VectoMatrix(Pos);
           VectorXm dE(3*NNode), dE_Damp(3*NNode);
           SparseMatrix<double> hEnergy1(3*NNode,3*NNode), hEnergy2(3*NNode,3*NNode),hDamp1(3*NNode,3*NNode), hDamp2(3*NNode,3*NNode);
	   vector<Matrix2m> atemp, btemp;
           dE_Damp.setZero();
           sm = new SimulationMesh(Vbar, V, F, t, constE, nu);
           if(DampingForce_Enabled && tcount >0){
                sm->testElasticEnergy(&dE_Damp, &hDamp1, &hDamp2, &atemp, &btemp, a_prev, b_prev);
           }
           sm->testElasticEnergy(&dE, &hEnergy1, &hEnergy2, &a_prev, &b_prev, abars, bbars);
           dE += cdamp1 / delt * dE_Damp;
	   Vel -= delt * MassInv * dE; 
           delete sm;
           //cout << "ForceImplicit" << endl << -0.5*dE << endl; 
           // cout << "ForceExplicit" << endl << Force() << endl; 
	}
        


	//----------------------Update with Verlocity Verlet---------------------------------------
	if(Implicit == 0){
	   V+=velocity*delt;
           MatrixXd FFtot;
           FFtot=Force();
//	   SimulationMesh *sm;
//           VectorXm dE(3*NNode);
//           SparseMatrix<double> hEnergy1(3*NNode,3*NNode), hEnergy2(3*NNode,3*NNode);
//           sm = new SimulationMesh(Vbar, V, F, t);
//           sm->testElasticEnergy(&dE, &hEnergy1, &hEnergy2, abars, bbars);
//	   if(tcount < 5){
//		cout << "Velet" << FFtot << endl << "Implicit Force" << endl << dE << endl;
//	   }
           
	   int j;
           for (j=0; j<NNode; j++){
             velocity.row(j)+=FFtot.row(j)*delt/Mass(j);
             if(gflag==1){ 
                velocity.col(1)-=gravity;
             }
           }
	}
	//------------------------------------------------------------------------------------------



	// Update the Moisture
	//Mtot_Mass=Mass_Total();
    	//Mtot_Stiffness=Stiffness_Total();
    	//Mtot_Stiffness*=DiffusionCoeff;
        if (tcount+InitialNum*storefreq >= (int) (sourcetime/delt)){
	    RHS_F.setZero();
	}
        //If a constant source exists, make the assumption that a constant moist is added to the vertices and diffuse afterwards.
        if(SourceType == 2){
            Moisture_v += Moist_source;
        }
        diffusion_prism();
        // Fix a few vertices from moving
        if (gflag==1){
 	    V.row(180)=V_ini1;
	    V.row(199)=V_ini2;
	}
}

MatrixXd VectoMatrix(VectorXd Vin){
	MatrixXd Mout(NNode,3);
	for (int i=0; i<3*NNode; i++){
	    Mout(i/3,i%3) = Vin(i);
	}
	return Mout;
}

VectorXd VMass(){
	int i;
        double areaT;
        Vector3d e1 , e2;
        VectorXd M(V.rows());
        M.setZero();
	for (i=0; i<F.rows(); i++){
        	e1 = V.row(F(i,1))-V.row(F(i,0));
        	e2 = V.row(F(i,2))-V.row(F(i,1));
		areaT=0.5*e1.cross(e2).norm();
                M(F(i,0))+=rho*areaT*t/3;
                M(F(i,1))+=rho*areaT*t/3;
                M(F(i,2))+=rho*areaT*t/3;
	}
	return M;
}

DiagonalMatrix<double,Dynamic> MasstoMatrix(){
        VectorXd MInv(3*NNode);
	for (int i=0; i<NNode; i++){
	    double temp = Mass(i);
	    MInv(3*i)=temp;
	    MInv(3*i+1)=temp;
	    MInv(3*i+2)=temp;
	}
	return MInv.asDiagonal();
}

DiagonalMatrix<double,Dynamic> InverseMass(){
        VectorXd MInv(3*NNode);
	for (int i=0; i<NNode; i++){
	    double temp = 1.0/Mass(i);
	    MInv(3*i)=temp;
	    MInv(3*i+1)=temp;
	    MInv(3*i+2)=temp;
	}
	return MInv.asDiagonal();
}

double TotEnergy( ){
	int i, j, k, NTri;
        double E;
	NTri=F.rows();
	MatrixXd FN(NTri,3), EN(3*NTri,3);
	E = 0;
	igl::per_face_normals(V, F, FN);
	EN=Enormal(FN);
	for (i=0; i<NTri; i++) {
		E+=Energy(constE, nu , t, i, V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)), EN.row(3*i), EN.row(3*i+1), EN.row(3*i+2));
	}
	return E;
}


//E1 computes the first energy term, v0 for rest 2D mesh, v for deformed 3D mesh, vbar for rest 3D vertices
double Energy(double constE, double nu, double t, int FaceIndex, Vector3d v1, Vector3d v2, Vector3d v3,Vector3d n1, Vector3d n2, Vector3d n3){
	Vector3d E;
	Matrix2d IA, IAbar, A, IB, IBbar, B;
	double dA, E1, E2, Etot;
	IA=I1(v1,v2,v3);
	IB=I2(v1,v2,v3,n1,n2,n3);
	IAbar=Ibartot.block(2*FaceIndex,0,2,2);
	IBbar=Ibartot.block(2*(FaceIndex+F.rows()),0,2,2);
	A=IAbar.inverse()*(IA-IAbar);
	B=IAbar.inverse()*(IB-IBbar);
	dA=0.5*sqrt(IAbar.determinant());
	E1=dA*t*(0.5*LameFirst*pow(A.trace(),2.0)+LameSecond*((A*A).trace()))/8;
	E2=dA*pow(t,3)*(0.5*LameFirst*pow(B.trace(),2.0)+LameSecond*((B*B).trace()))/24;	
	Etot=E1+E2;
        E << E1, E2, Etot;
	return Etot; 
}

// Return the force on every vertices
MatrixXd Force(){
	RowVector3d dval1, dval2, dval3;
	Matrix2d IA, IAbar, A, IB, IBbar, B, Rot, tmp, Inv;
	Matrix3d Ed;
	int i, NTri;
	double c1, c2, dval, dA;
	NTri=F.rows();
	MatrixXd dN, dI(18,2), Tr(2,3), TrDamp(2,3), FF(NNode,3), FF2(NNode,3), FF1(NNode,3), FN(NTri,3), EN(3*NTri,3), FDamp1(NNode,3), FDamp2(NNode,3);
	igl::per_face_normals(V, F, FN);
	EN=Enormal(FN);	
	dN=dEnormM(); //dN is the derivative of egde normal w.r.t. the four vectices cooresponding to the edge
        FF2.setZero();
	FF1.setZero();
	FDamp1.setZero();
	FDamp2.setZero();
	for (i=0; i<NTri; i++) {
		Ed << V.row(F(i,1))-V.row(F(i,0)), V.row(F(i,2))-V.row(F(i,1)), V.row(F(i,0))-V.row(F(i,2)); //Ed.row(i) is the ith  edge vector for deformed triangle
		IB=I2(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),EN.row(3*i),EN.row(3*i+1),EN.row(3*i+2));
		IA=I1(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
		IAbar=Ibartot.block(2*i,0,2,2);
		IBbar=Ibartot.block(2*(i+F.rows()),0,2,2);
		Inv=IAbar.inverse();
                B=Inv*(IB-IBbar);
		A=Inv*(IA-IAbar);
		/*
		if(tcount<2){
		  cout << "Face " << i << endl;
		  cout <<  "IAbar" << endl << IAbar << endl << "IA" << endl << IA << endl;
		  cout <<  "IBbar" << endl << IBbar << endl << "IB" << endl << IB << endl;
		}
		*/
		dI=DelI(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
		dA=0.5*sqrt(IAbar.determinant());
		c1=0.5*LameFirst*dA/24;
		c2=LameSecond*dA/24;
		Tr << B.trace()*(Inv*CanonM1).trace(), B.trace()*(Inv*CanonM2).trace(), B.trace()*(Inv*CanonM3).trace(), (B*Inv*CanonM1).trace(), (B*Inv*CanonM2).trace(), (B*Inv*CanonM3).trace();
		Matrix2d BDamp;
		BDamp = Inv*(IB-IBold.block(2*i,0,2,2));
		TrDamp << BDamp.trace()*(Inv*CanonM1).trace(), BDamp.trace()*(Inv*CanonM2).trace(), BDamp.trace()*(Inv*CanonM3).trace(), (BDamp*Inv*CanonM1).trace(), (BDamp*Inv*CanonM2).trace(), (BDamp*Inv*CanonM3).trace();
		for (int j=0; j<3; j++){
			for (int k=0; k<4; k++){
				if (EdgeNV(3*i+j,k)!=-1){
					dval1=2*Ed.row((j+1)%3)*dN.block(9*i+3*j,3*k,3,3);	// orginate from the term 2<dn(j),e(j+1)>, dn(j)(dV(EdgeNV(j,k)))
					dval2=-2*Ed.row((j+2)%3)*dN.block(9*i+3*j,3*k,3,3); // originate form the term -2<dn(j),e(j-1)>
					FF2.row(EdgeNV(3*i+j,k))+= -2*pow(t,3)*(c1*(Tr(0,(j+1)%3)-Tr(0,j)-Tr(0,(j+2)%3))+c2*(Tr(1,(j+1)%3)-Tr(1,j)-Tr(1,(j+2)%3)))*dval1;
					FF2.row(EdgeNV(3*i+j,k))+= -2*pow(t,3)*(c1*(Tr(0,(j+2)%3)-Tr(0,j)-Tr(0,(j+1)%3))+c2*(Tr(1,(j+2)%3)-Tr(1,j)-Tr(1,(j+1)%3)))*dval2;
               				if(DampingForce_Enabled && tcount >0){
						FDamp2.row(EdgeNV(3*i+j,k))+= -2*pow(t,3)*(c1*(TrDamp(0,(j+1)%3)-TrDamp(0,j)-TrDamp(0,(j+2)%3))+c2*(TrDamp(1,(j+1)%3)-TrDamp(1,j)-TrDamp(1,(j+2)%3)))*dval1;
						FDamp2.row(EdgeNV(3*i+j,k))+= -2*pow(t,3)*(c1*(TrDamp(0,(j+2)%3)-TrDamp(0,j)-TrDamp(0,(j+1)%3))+c2*(TrDamp(1,(j+2)%3)-TrDamp(1,j)-TrDamp(1,(j+1)%3)))*dval2;
					}		
				}
			}
			dval3=2*(EN.row(3*i+((j+2)%3))-EN.row(3*i+((j+1)%3))); // <n3-n2,dv1>, <n3-n3,dv2>
			FF2.row(F(i,j))+= 2*pow(t,3)*(c1*(Tr(0,j)-Tr(0,(j+1)%3)-Tr(0,(j+2)%3))+c2*(Tr(1,j)-Tr(1,(j+1)%3)-Tr(1,(j+2)%3)))*dval3;// dQ(j)(dV(j))
			FF2.row(F(i,(j+1)%3))+= -2*pow(t,3)*(c1*(Tr(0,j)-Tr(0,(j+1)%3)-Tr(0,(j+2)%3))+c2*(Tr(1,j)-Tr(1,(j+1)%3)-Tr(1,(j+2)%3)))*dval3;//dQ(j)(dV(j+1))
               		if(DampingForce_Enabled && tcount >0){
				FDamp2.row(F(i,j))+= 2*pow(t,3)*(c1*(TrDamp(0,j)-TrDamp(0,(j+1)%3)-TrDamp(0,(j+2)%3))+c2*(TrDamp(1,j)-TrDamp(1,(j+1)%3)-TrDamp(1,(j+2)%3)))*dval3;// dQ(j)(dV(j))
				FDamp2.row(F(i,(j+1)%3))+= -2*pow(t,3)*(c1*(TrDamp(0,j)-TrDamp(0,(j+1)%3)-TrDamp(0,(j+2)%3))+c2*(TrDamp(1,j)-TrDamp(1,(j+1)%3)-TrDamp(1,(j+2)%3)))*dval3;//dQ(j)(dV(j+1))
			}
		}
		Matrix2d ADamp;
		ADamp = Inv*(IA-IAold.block(2*i,0,2,2));
		for (int j=0; j<3; j++){
			for(int k=0; k<3; k++){
				tmp=Inv*dI.block(6*j+2*k,0,2,2);
				FF1(F(i,j),k)+=-6*t*(c1*A.trace()*tmp.trace()+c2*(A*tmp).trace()); // Force due to the fisrt fundamental term
               			if(DampingForce_Enabled && tcount >0){
					FDamp1(F(i,j),k)+=-6*t*(c1*ADamp.trace()*tmp.trace()+c2*(ADamp*tmp).trace()); // Damping force for the streching energy
				}
                /*if((A*tmp).trace()>0.000001){
		cout << "Vertex " << F(i,j) << " " << k << " direction" << endl << "tcount " << tcount << endl;
	        cout << "A " << endl << A << endl << "Inv" << endl << Inv << endl;
         	cout << "IA-IAbar_" << endl << IA-IAbar << endl;
         	}
		*/
                
			}	
		}
		/*
               //Calculate Damping Force 
               if(DampingForce_Enabled && tcount >0){
		Matrix2d ADamp, BDamp;
		ADamp = Inv*(IA-IAold.block(2*i,0,2,2));
		BDamp = Inv*(IB-IBold.block(2*i,0,2,2));
		for (int j=0; j<3; j++){
			for(int k=0; k<3; k++){
				// Fdamp1 is for Damping force due to stretching
				tmp=Inv*dI.block(6*j+2*k,0,2,2);
				FDamp1(F(i,j),k)+=-6*t*(c1*ADamp.trace()*tmp.trace()+c2*(ADamp*tmp).trace()); // Damping force for the streching energy
				// Fdamp1 is for Damping force due to stretching
			}	
		}
		Tr << BDamp.trace()*(Inv*CanonM1).trace(), BDamp.trace()*(Inv*CanonM2).trace(), BDamp.trace()*(Inv*CanonM3).trace(), (BDamp*Inv*CanonM1).trace(), (BDamp*Inv*CanonM2).trace(), (BDamp*Inv*CanonM3).trace();
		for (int j=0; j<3; j++){
			for (int k=0; k<4; k++){
				if (EdgeNV(3*i+j,k)!=-1){
					dval1=2*Ed.row((j+1)%3)*dN.block(9*i+3*j,3*k,3,3);	// orginate from the term 2<dn(j),e(j+1)>, dn(j)(dV(EdgeNV(j,k)))
					dval2=-2*Ed.row((j+2)%3)*dN.block(9*i+3*j,3*k,3,3); // originate form the term -2<dn(j),e(j-1)>
					FDamp2.row(EdgeNV(3*i+j,k))+= -2*pow(t,3)*(c1*(Tr(0,(j+1)%3)-Tr(0,j)-Tr(0,(j+2)%3))+c2*(Tr(1,(j+1)%3)-Tr(1,j)-Tr(1,(j+2)%3)))*dval1;
					FDamp2.row(EdgeNV(3*i+j,k))+= -2*pow(t,3)*(c1*(Tr(0,(j+2)%3)-Tr(0,j)-Tr(0,(j+1)%3))+c2*(Tr(1,(j+2)%3)-Tr(1,j)-Tr(1,(j+1)%3)))*dval2;
				}
			}
			dval3=2*(EN.row(3*i+((j+2)%3))-EN.row(3*i+((j+1)%3))); // <n3-n2,dv1>, <n3-n3,dv2>
			FDamp2.row(F(i,j))+= 2*pow(t,3)*(c1*(Tr(0,j)-Tr(0,(j+1)%3)-Tr(0,(j+2)%3))+c2*(Tr(1,j)-Tr(1,(j+1)%3)-Tr(1,(j+2)%3)))*dval3;// dQ(j)(dV(j))
			FDamp2.row(F(i,(j+1)%3))+= -2*pow(t,3)*(c1*(Tr(0,j)-Tr(0,(j+1)%3)-Tr(0,(j+2)%3))+c2*(Tr(1,j)-Tr(1,(j+1)%3)-Tr(1,(j+2)%3)))*dval3;//dQ(j)(dV(j+1))
		}
		// FDamp should be divived by delt since it is related to the strain rate, however we absorb the term into cdamp
	       }
		*/
		IAold.block(2*i,0,2,2)=IA;
		IBold.block(2*i,0,2,2)=IB;
	}
	double C;
	C=-1/(8*area*area);
	FF2=-C*FF2;
	FDamp1 *= cdamp1/delt;
	FDamp2 *= cdamp2*-C/delt;
	FF = FF1+FF2+FDamp1+FDamp2;
	/*// For Testing
	if(tcount < 2){
        cout << tcount << endl << "FF1" << endl << FF1 << endl;
        cout << tcount << endl << "FF2" << endl << FF2 << endl;
	}
        */
	return FF;
}

MatrixXd DelI(Vector3d v1, Vector3d v2, Vector3d v3){ // return delta I w.r.t to the change in Vi, V2, V3
	MatrixXd T1(6,2), T2(6,2), T3(6,2), dI(18,2);
	double c;
	int i;
	for (i=0; i<3; i++){
		T1.block(2*i,0,2,2)=2*(v3(i)-v2(i))*CanonM1+2*(v3(i)+v2(i)-2*v1(i))*CanonM2+2*(v2(i)-v3(i))*CanonM3;
		T2.block(2*i,0,2,2)=2*(v3(i)-v1(i))*CanonM1+2*(v1(i)-v3(i))*CanonM2+2*(v1(i)+v3(i)-2*v2(i))*CanonM3;
		T3.block(2*i,0,2,2)=2*(v1(i)+v2(i)-2*v3(i))*CanonM1+2*(v1(i)-v2(i))*CanonM2+2*(v2(i)-v1(i))*CanonM3;
	}
	dI << T1, T2, T3;
	c=-1/(8*area*area);
	dI=c*dI;
	return dI;
}

Vector2d FoldCircle(double mid, double xpos, double r){
	Vector2d vec;
	double theta;
	theta=(xpos-mid)/r;
	vec << mid+r*sin(theta), r*cos(theta);
	return vec;
}

VectorXi NeighborF( ){
	int NTri=F.rows(),i,j,m,n,flag;
	VectorXi EdgeFf(3*NTri);	// EdgeF records the neighboring face of an edge in the tri. mesh in ccw order, #n edge in  face #i is Edge(3*i+n) 
	EdgeFf=-1*EdgeFf.setOnes(3*NTri);
	for (i=0; i< NTri; i++){
		for (j=0; j<3; j++){
			m=0;
			flag=0;
			while (m<NTri && flag<2){
				flag=0;
				if (m!=i){	//avoid self counting
					for(n=0; n<3; n++){
						if (F(i,j)==F(m,n) || F(i,(j+1)%3)==F(m,n)){
							flag++;
						}
					}
				}
				m++;
			}
		
			if (flag==2){
				EdgeFf(3*i+j)=m-1;
			}
		}
	}
	return EdgeFf;
}

MatrixXi VofEdgeN( ){  // the related vertices to an edge normal
	int NEdge=EdgeF.size(), i ,j, NTri=F.rows(), m;
	MatrixXi EdgeNVf(NEdge,4);
	EdgeNVf.col(3).setConstant(-1);
	for (i=0; i<NTri; i++){
		for (j=0; j<3; j++){
			EdgeNVf(3*i+j,0)= F(i,j);
			EdgeNVf(3*i+j,1)= F(i,(j+1)%3);
			EdgeNVf(3*i+j,2)= F(i,(j+2)%3);
			if (EdgeF(3*i+j)!= -1){
				for (m=0; m<3; m++){
					if(F(i,j)!=F(EdgeF(3*i+j),m)&&F(i,(j+1)%3)!=F(EdgeF(3*i+j),m)){
						EdgeNVf(3*i+j,3)=F(EdgeF(3*i+j),m);
					}
				}
			}	
		}
	} 
	return EdgeNVf;
}

Matrix3d CrossM(Vector3d v){ // return the matrix A such that axb=Ab, A is the cross product matrix
	Matrix3d M;
	M << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
	return M;
}

Matrix3d dNormVecM(Vector3d v){ // return the matrix M such that for F(x)=x/|x|, dF=M*dx
	Matrix3d M, I;
	double a=1/v.norm();
	I.setIdentity();
	M=a*I-pow(a,3)*v*v.transpose();
	return M;
}

MatrixXd dEnormM( ){ // return the derivative of each edge normal w. the vertices related to the normal of edge
	// EdgeNV is the matrix that contains the related vertices to an edge, matrix V is just the vertice coordinate matrix
	int i, NEdge=EdgeNV.rows();
	Matrix3d Mn, M1, M2, c1_1, c2_1, c1_2, c2_2;
	MatrixXd M(3*NEdge,12);   // M.block(3*i,3*j,3,3) corresponse to the derivative of the i-th edge with EdgeV(i,j) component
	Vector3d e1_1, e2_1, e1_2, e2_2;
	for (i=0; i<NEdge; i++){
		e1_1=V.row(EdgeNV(i,1))-V.row(EdgeNV(i,0));   // edg1 1 in face 1 (current triangle)
		e2_1=V.row(EdgeNV(i,2))-V.row(EdgeNV(i,1));
		e1_2=V.row(EdgeNV(i,0))-V.row(EdgeNV(i,1));   // edge 1 in face 2 (neighboring triangle)
		if (EdgeNV(i,3)==-1){	// no neighboring face for Edge(i)
			Mn.setIdentity();
			M2.setZero();
			e2_2.setZero();
		}
		else {
			e2_2=V.row(EdgeNV(i,3))-V.row(EdgeNV(i,0));
			Mn=dNormVecM(e1_1.cross(e2_1).normalized()+e1_2.cross(e2_2).normalized());
			M2=dNormVecM(e1_2.cross(e2_2));
		}
		c1_1=CrossM(e1_1);
		c2_1=CrossM(e2_1);
		c1_2=CrossM(e1_2);
		c2_2=CrossM(e2_2);
		M1=dNormVecM(e1_1.cross(e2_1));
		M.block(3*i,0,3,3)=Mn*(M1*c2_1-M2*(c1_2+c2_2));
		M.block(3*i,3,3,3)=Mn*(-M1*(c1_1+c2_1)+M2*c2_2);
		M.block(3*i,6,3,3)=Mn*M1*c1_1;
		M.block(3*i,9,3,3)=Mn*M2*c1_2;
	}
	return M;
}
//Enormal calculates the normal of the edges by averaging over neighboring faces and the #n edge in #i face has normal EdgeN(3*i+n)
MatrixXd Enormal(MatrixXd FNf){	// FN contains the face normals and EdgeF contains the neighboring faces of each edge
	MatrixXd EdgeN(EdgeF.size(),3);
	int i;
	for(i=0; i<EdgeF.size(); i++){
		if (EdgeF(i)==-1){
			EdgeN.row(i) << FNf.row(i/3);
		}
		else {
			EdgeN.row(i)=FNf.row(i/3)+FNf.row(EdgeF(i));
			EdgeN.row(i)=EdgeN.row(i).normalized();
		}	
	}
	return EdgeN;
} 

//I1 computes the matrix of the first fundamental form based on the vertices of the triangle mesh
Matrix2d I1(Vector3d v1d, Vector3d v2d, Vector3d v3d){
	double q1, q2, q3;
	Matrix2d I;
	q1=(v2d-v1d).squaredNorm();
	q2=(v3d-v2d).squaredNorm();
	q3=(v1d-v3d).squaredNorm();
	I=QMatrix(q1,q2,q3);
	return I;	
}

Matrix2d I2(Vector3d v1, Vector3d v2, Vector3d v3, Vector3d n1, Vector3d n2, Vector3d n3){
	Vector3d e1, e2, e3;
	double q1, q2, q3;
	Matrix2d I;
	e1=v2-v1;
	e2=v3-v2;
	e3=v1-v3;
	q1=2*e1.dot(n2-n3);
	q2=2*e2.dot(n3-n1);
	q3=2*e3.dot(n1-n2);
	I=QMatrix(q1,q2,q3);
	return I;
}

//QMatrix returns the quadratic function Q in discrete triangular mesh as a matrix (operator)
Matrix2d QMatrix(double q1, double q2, double q3){
	Matrix2d Q;
	double c;
	c=-1/(8*area*area);
	Q=c*((q1-q2-q3)*CanonM1+(q2-q3-q1)*CanonM2+(q3-q1-q2)*CanonM3);
	return Q;
}


//edge computes the edge vector given the vertices of the triangle   
Vector3d edge(Vector3d a, Vector3d b){
        Vector3d c=b-a;
        return c;
}

//rMatrix computes the rotation matrix to generate the dual edges, e1 e2 in ccw ordering
Matrix3d rMatrix(Vector3d e1, Vector3d e2){
	Vector3d axis=e1.cross(e2);
	//axis=axis.normalize();
	Matrix3d rot;	//rot is the planar rotation matrix wrt to the axis
	rot=AngleAxisd(-0.5*M_PI,axis.normalized());
	return rot;
}

//Diffusion in 3D using a prism element
void diffusion_prism(){
	//SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        SimplicialLDLT<SparseMatrix<double> >  solver;
	solver.compute(Mtot_Mass/delt+Mtot_Stiffness);
	if(solver.info()!=Success) {
  		// decomposition failed
  	   return;
        }
	Moisture_v = solver.solve(Mtot_Mass*Moisture_v/delt+RHS_F);	
}

// Calculate the Jacobian matrix of Face i
Matrix3d M_Jacobian(int i){
	Matrix3d J;
	Vector3d v1,v2,v3,n;
	v1=Vbar.row(F(i,0));
	v2=Vbar.row(F(i,1));
	v3=Vbar.row(F(i,2));
	n = ((v2-v1).cross(v3-v1)).normalized();
	J << -v1(0)+v2(0), -v1(0)+v3(0), 0, -v1(1)+v2(1), -v1(1)+v3(1), 0, -v1(2)+v2(2), -v1(2)+v3(2), 0;
	J.col(2)= t*n/2;
	return J; 
}

//Calculate the RHS of the diffusion eq. with constant rate of dissipation which is linear from ConstS on the top surface
//to zero in the bottom surface
Vector6d UniTopF_prism(int i){
	Vector6d RHS_v;
	RHS_v << 2.0/9, 2.0/9, 2.0/9, 1.0/9, 1.0/9, 1.0/9;
	return Dissipation*M_Jacobian(i).determinant()*RHS_v/2;
}

//TODO:Wrong!!! Have to recalculate the value based on the integration in the prism 
//o.w. neighboring vertices with 0 contribution will lead to negative moisture change
//Appling a plane source from the bottom surface. i.e. soaking in water
/*
Vector6d PlaneSourceBot(int i){
	Vector6d RHS_v;
        RHS_v.setZero();
        for(int j=0; j<3; j++){
	   // Assuming that V is the coordinate of the mid surface
           if(V(F(i,j),2)-t/2<SourceZ+1e-9){
	      	RHS_v(j)=1.0/9;
	      	RHS_v(j+3)=2.0/9;
	   }
	}
	return Dissipation*M_Jacobian(i).determinant()*RHS_v/2;
}
*/

VectorXd F_Total(){
	VectorXd Tot_F(2*NNode);
	Vector6d Vec_F;
	Tot_F.setZero();
	int NTri, index1, index2, index3, index4, index5, index6;
	NTri=F.rows();
	if(SourceType != 0 && SourceType != 2){
	   for(int i=0; i<NTri; i++){
	   	Vec_F=funcRHS[SourceType-1](i);
	   	index1=F(i,0);
	   	index2=F(i,1);
	   	index3=F(i,2);
	   	index4=index1+NNode;// Should be equalvent to index1+NNode
	   	index5=index2+NNode;
	   	index6=index3+NNode;
	   	Tot_F(index1)+=Vec_F(0);
	   	Tot_F(index2)+=Vec_F(1);
	   	Tot_F(index3)+=Vec_F(2);
	   	Tot_F(index4)+=Vec_F(3);
	   	Tot_F(index5)+=Vec_F(4);
	   	Tot_F(index6)+=Vec_F(5);
	   }
	}
	return Tot_F; 
}

// Calculate the Mass matrix for prism i
Matrix6d Mass_prism(int i){
	Matrix6d M_Mass;
	M_Mass << 1.0/18, 1.0/36, 1.0/36, 1.0/36, 1.0/72, 1.0/72, 1.0/36, 1.0/18, 1.0/36, 1.0/72, 1.0/36, 1.0/72, 1.0/36, 1.0/36, 1.0/18, 1.0/72, 1.0/72, 1.0/36, 1.0/36, 1.0/72, 1.0/72, 1.0/18, 1.0/36, 1.0/36, 1.0/72, 1.0/36, 1.0/72, 1.0/36, 1.0/18, 1.0/36, 1.0/72, 1.0/72, 1.0/36, 1.0/36, 1.0/36, 1.0/18; 
	return M_Jacobian(i).determinant()*M_Mass; 	
}

SparseMatrix<double> Mass_Total(){
	SparseMatrix<double> Tot_Mass(2*NNode,2*NNode);
	Matrix6d M_Mass;
	Tot_Mass.setZero();
	int NTri, index1, index2, index3, index4, index5, index6;
	NTri=F.rows();
	for(int i=0; i<NTri;i++){
		M_Mass=Mass_prism(i);
		index1=F(i,0);
		index2=F(i,1);
		index3=F(i,2);
		index4=index1+NNode;// Should be equalvent to index1+NNode
		index5=index2+NNode;
		index6=index3+NNode;
		Tot_Mass.coeffRef(index1,index1)+=M_Mass(0,0);
		Tot_Mass.coeffRef(index1,index2)+=M_Mass(0,1);
		Tot_Mass.coeffRef(index1,index3)+=M_Mass(0,2);
		Tot_Mass.coeffRef(index1,index4)+=M_Mass(0,3);
		Tot_Mass.coeffRef(index1,index5)+=M_Mass(0,4);
		Tot_Mass.coeffRef(index1,index6)+=M_Mass(0,5);
		Tot_Mass.coeffRef(index2,index1)+=M_Mass(1,0);
		Tot_Mass.coeffRef(index2,index2)+=M_Mass(1,1);
		Tot_Mass.coeffRef(index2,index3)+=M_Mass(1,2);
		Tot_Mass.coeffRef(index2,index4)+=M_Mass(1,3);
		Tot_Mass.coeffRef(index2,index5)+=M_Mass(1,4);
		Tot_Mass.coeffRef(index2,index6)+=M_Mass(1,5);
		Tot_Mass.coeffRef(index3,index1)+=M_Mass(2,0);
		Tot_Mass.coeffRef(index3,index2)+=M_Mass(2,1);
		Tot_Mass.coeffRef(index3,index3)+=M_Mass(2,2);
		Tot_Mass.coeffRef(index3,index4)+=M_Mass(2,3);
		Tot_Mass.coeffRef(index3,index5)+=M_Mass(2,4);
		Tot_Mass.coeffRef(index3,index6)+=M_Mass(2,5);
		Tot_Mass.coeffRef(index4,index1)+=M_Mass(3,0);
		Tot_Mass.coeffRef(index4,index2)+=M_Mass(3,1);
		Tot_Mass.coeffRef(index4,index3)+=M_Mass(3,2);
		Tot_Mass.coeffRef(index4,index4)+=M_Mass(3,3);
		Tot_Mass.coeffRef(index4,index5)+=M_Mass(3,4);
		Tot_Mass.coeffRef(index4,index6)+=M_Mass(3,5);
		Tot_Mass.coeffRef(index5,index1)+=M_Mass(4,0);
		Tot_Mass.coeffRef(index5,index2)+=M_Mass(4,1);
		Tot_Mass.coeffRef(index5,index3)+=M_Mass(4,2);
		Tot_Mass.coeffRef(index5,index4)+=M_Mass(4,3);
		Tot_Mass.coeffRef(index5,index5)+=M_Mass(4,4);
		Tot_Mass.coeffRef(index5,index6)+=M_Mass(4,5);
		Tot_Mass.coeffRef(index6,index1)+=M_Mass(5,0);
		Tot_Mass.coeffRef(index6,index2)+=M_Mass(5,1);
		Tot_Mass.coeffRef(index6,index3)+=M_Mass(5,2);
		Tot_Mass.coeffRef(index6,index4)+=M_Mass(5,3);
		Tot_Mass.coeffRef(index6,index5)+=M_Mass(5,4);
		Tot_Mass.coeffRef(index6,index6)+=M_Mass(5,5);
	}
	return Tot_Mass;
}                                                
                                                 
// Calculate the stiffness matrix for prism i   
Matrix6d Stiffness_prism(int i){                
	Matrix6d M_Stiffness;                   
        Matrix3d J,J_invT;
	J=M_Jacobian(i);
	J_invT=(J.inverse()).transpose();
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	double K11, K12, K13, K14, K15, K16, K22, K23, K24, K25, K26, K33, K34, K35, K36, K44, K45, K46, K55, K56, K66;
	a11=J_invT(0,0);
	a12=J_invT(0,1);
	a13=J_invT(0,2);
	a21=J_invT(1,0);
	a22=J_invT(1,1);
	a23=J_invT(1,2);
	a31=J_invT(2,0);
	a32=J_invT(2,1);
	a33=J_invT(2,2);
        K11= (8*a11*a11+16*a11*a12+8*a12*a12-4*a11*a13-4*a12*a13+a13*a13+8*a21*a21+16*a21*a22+8*a22*a22-4*a21*a23-4*a22*a23+a23*a23+8*a31*a31+16*a31*a32+8*a32*a32-4*a31*a33-4*a32*a33+a33*a33)/24;
        K12= (-16*a11*a11-16*a11*a12-4*a12*a13+a13*a13-16*a21*a21-16*a21*a22-4*a22*a23+a23*a23-16*a31*a31-16*a31*a32-4*a32*a33+a33*a33)/48;
        K13=1.0/48*(-16*a12*a12 +a13*a13 -4*a11*(4*a12+a13)-16*a21*a22-16*a22*a22-4*a21*a23+a23*a23-16*a31*a32-16*a32*a32-4*a31*a33+a33*a33);
        K14= (4*a11*a11+8*a11*a12+4*a12*a12-a13*a13+4*a21*a21+8*a21*a22+4*a22*a22-a23*a23+4*a31*a31+8*a31*a32+4*a32*a32-a33*a33)/24;
	K15=1.0/48*(-8*a11*a11 - 8*a11*(a12 - a13) + 4*a12*a13 - a13*a13 - 8*a21*a21 -8*a21*a22+8*a21*a23+4*a22*a23-a23*a23-8*a31*a31-8*a31*a32+8*a31*a33+4*a32*a33-a33*a33);
	K16=1.0/48*(-8*a11*a12-8*a12*a12+4*a11*a13+8*a12*a13-a13*a13-8*a21*a22-8*a22*a22+4*a21*a23+8*a22*a23-a23*a23-8*a31*a32-8*a32*a32+4*a31*a33+8*a32*a33-a33*a33);
	K22=1.0/24*(8*a11*a11+4*a11*a13+a13*a13+8*a21*a21+4*a21*a23+a23*a23+8*a31*a31+4*a31*a33+a33*a33);
	K23=1.0/48*(4*a12*a13+a13*a13+4*a11*(4*a12+a13)+16*a21*a22+4*a21*a23+4*a22*a23+a23*a23+16*a31*a32+4*a31*a33+4*a32*a33+a33*a33);
	K24=1.0/48*(-8*a11*a11-4*a12*a13-a13*a13-8*a11*(a12+a13)-8*a21*a21-8*a21*a22-8*a21*a23-4*a22*a23-a23*a23-8*a31*a31-8*a31*a32-8*a31*a33-4*a32*a33-a33*a33);
	K25=1.0/24*(4*a11*a11-a13*a13+4*a21*a21-a23*a23+4*a31*a31-a33*a33);
	K26=1.0/48*(8*a11*a12-4*a11*a13+4*a12*a13-a13*a13+8*a21*a22-4*a21*a23+4*a22*a23-a23*a23+8*a31*a32-4*a31*a33+4*a32*a33-a33*a33);
	K33=1.0/24*(8*a12*a12+4*a12*a13+a13*a13+8*a22*a22+4*a22*a23+a23*a23+8*a32*a32+4*a32*a33+a33*a33);
	K34=1.0/48*(-8*a12*a12-8*a12*a13-a13*a13-4*a11*(2*a12+a13)-8*a21*a22-8*a22*a22-4*a21*a23-8*a22*a23-a23*a23-8*a31*a32-8*a32*a32-4*a31*a33-8*a32*a33-a33*a33);
	K35=1.0/48*(-4*a12*a13-a13*a13+4*a11*(2*a12+a13)+8*a21*a22+4*a21*a23-4*a22*a23-a23*a23+8*a31*a32+4*a31*a33-4*a32*a33-a33*a33);
	K36=-1.0/24*(-4*a12*a12+a13*a13-4*a22*a22+a23*a23-4*a32*a32+a33*a33);
	K44=1.0/24*(8*a11*a11+8*a12*a12+4*a12*a13+a13*a13+4*a11*(4*a12+a13)+8*a21*a21+16*a21*a22+8*a22*a22+4*a21*a23+4*a22*a23+a23*a23+8*a31*a31+16*a31*a32+8*a32*a32+4*a31*a33+4*a32*a33+a33*a33);
	K45=1.0/48*(-16*a11*a11-16*a11*a12+4*a12*a13+a13*a13-16*a21*a21-16*a21*a22+4*a22*a23+a23*a23-16*a31*a31-16*a31*a32+4*a32*a33+a33*a33);
	K46=1.0/48*(-16*a12*a12+a13*a13+4*a11*(-4*a12+a13)-16*a21*a22-16*a22*a22+4*a21*a23+a23*a23-16*a31*a32-16*a32*a32+4*a31*a33+a33*a33);
	K55=1.0/24*(8*a11*a11-4*a11*a13+a13*a13+8*a21*a21-4*a21*a23+a23*a23+8*a31*a31-4*a31*a33+a33*a33);
	K56=1.0/48*(16*a11*a12-4*a11*a13-4*a12*a13+a13*a13+16*a21*a22-4*a21*a23-4*a22*a23+a23*a23+16*a31*a32-4*a31*a33-4*a32*a33+a33*a33);
	K66=1.0/24*(8*a12*a12-4*a12*a13+a13*a13+8*a22*a22-4*a22*a23+a23*a23+8*a32*a32-4*a32*a33+a33*a33);
	M_Stiffness << K11, K12, K13, K14, K15, K16, K12, K22, K23, K24, K25, K26, K13, K23, K33, K34, K35, K36, K14, K24, K34, K44, K45, K46, K15, K25, K35, K45, K55, K56, K16, K26, K36, K46, K56, K66;
	return J.determinant()*M_Stiffness;
}

SparseMatrix<double> Stiffness_Total(){
	SparseMatrix<double> Tot_Stiffness(2*NNode,2*NNode);
	Matrix6d M_Stiff;
	Tot_Stiffness.setZero();
	int NTri, index1, index2, index3, index4, index5, index6;
	NTri=F.rows();
	for(int i=0; i<NTri;i++){
		M_Stiff=Stiffness_prism(i);
		index1=F(i,0);
		index2=F(i,1);
		index3=F(i,2);
		index4=index1+NNode;// Should be equalvent to index1+NNode
		index5=index2+NNode;
		index6=index3+NNode;
		Tot_Stiffness.coeffRef(index1,index1)+=M_Stiff(0,0);
		Tot_Stiffness.coeffRef(index1,index2)+=M_Stiff(0,1);
		Tot_Stiffness.coeffRef(index1,index3)+=M_Stiff(0,2);
		Tot_Stiffness.coeffRef(index1,index4)+=M_Stiff(0,3);
		Tot_Stiffness.coeffRef(index1,index5)+=M_Stiff(0,4);
		Tot_Stiffness.coeffRef(index1,index6)+=M_Stiff(0,5);
		Tot_Stiffness.coeffRef(index2,index1)+=M_Stiff(1,0);
		Tot_Stiffness.coeffRef(index2,index2)+=M_Stiff(1,1);
		Tot_Stiffness.coeffRef(index2,index3)+=M_Stiff(1,2);
		Tot_Stiffness.coeffRef(index2,index4)+=M_Stiff(1,3);
		Tot_Stiffness.coeffRef(index2,index5)+=M_Stiff(1,4);
		Tot_Stiffness.coeffRef(index2,index6)+=M_Stiff(1,5);
		Tot_Stiffness.coeffRef(index3,index1)+=M_Stiff(2,0);
		Tot_Stiffness.coeffRef(index3,index2)+=M_Stiff(2,1);
		Tot_Stiffness.coeffRef(index3,index3)+=M_Stiff(2,2);
		Tot_Stiffness.coeffRef(index3,index4)+=M_Stiff(2,3);
		Tot_Stiffness.coeffRef(index3,index5)+=M_Stiff(2,4);
		Tot_Stiffness.coeffRef(index3,index6)+=M_Stiff(2,5);
		Tot_Stiffness.coeffRef(index4,index1)+=M_Stiff(3,0);
		Tot_Stiffness.coeffRef(index4,index2)+=M_Stiff(3,1);
		Tot_Stiffness.coeffRef(index4,index3)+=M_Stiff(3,2);
		Tot_Stiffness.coeffRef(index4,index4)+=M_Stiff(3,3);
		Tot_Stiffness.coeffRef(index4,index5)+=M_Stiff(3,4);
		Tot_Stiffness.coeffRef(index4,index6)+=M_Stiff(3,5);
		Tot_Stiffness.coeffRef(index5,index1)+=M_Stiff(4,0);
		Tot_Stiffness.coeffRef(index5,index2)+=M_Stiff(4,1);
		Tot_Stiffness.coeffRef(index5,index3)+=M_Stiff(4,2);
		Tot_Stiffness.coeffRef(index5,index4)+=M_Stiff(4,3);
		Tot_Stiffness.coeffRef(index5,index5)+=M_Stiff(4,4);
		Tot_Stiffness.coeffRef(index5,index6)+=M_Stiff(4,5);
		Tot_Stiffness.coeffRef(index6,index1)+=M_Stiff(5,0);
		Tot_Stiffness.coeffRef(index6,index2)+=M_Stiff(5,1);
		Tot_Stiffness.coeffRef(index6,index3)+=M_Stiff(5,2);
		Tot_Stiffness.coeffRef(index6,index4)+=M_Stiff(5,3);
		Tot_Stiffness.coeffRef(index6,index5)+=M_Stiff(5,4);
		Tot_Stiffness.coeffRef(index6,index6)+=M_Stiff(5,5);
	}
	return Tot_Stiffness;
}                                                
