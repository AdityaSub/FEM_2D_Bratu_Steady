#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
#include<math.h>
#include<cmath>
#include "Mesh.h"

#define PI 3.14159265358979323846

using namespace std;
using namespace Eigen;

extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);

// constructor
Mesh::Mesh(string& fileName) : meshFileName(fileName) {
	readMesh();
	for (size_t i = 0; i < static_cast<size_t>(elements.size()); i++) {
		const Node n1(elements[i][0], nodes[elements[i][0]][0], nodes[elements[i][0]][1]);
		const Node n2(elements[i][1], nodes[elements[i][1]][0], nodes[elements[i][1]][1]);
		const Node n3(elements[i][2], nodes[elements[i][2]][0], nodes[elements[i][2]][1]);
		mesh.push_back(Element(n1, n2, n3, i));		
	}
	cout << "Grid generated! Bounds: x_min = " << x_min << ", y_min = " << y_min << ", x_max = " << x_max << ", y_max = " << y_max << ", nodes = " << nodes.size() << ", elements = " << elements.size() << endl;
	
	// calculate element stiffnesses (to be used in Jacobian and residual calculations later)
	for (vector<Element>::iterator it = mesh.begin(); it != mesh.end(); it++) 
		it->calcElementStiffness();  
	
	MatCreate(PETSC_COMM_WORLD,&J);
    MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,nodes.size(),nodes.size());
    MatSetUp(J);
    MatZeroEntries(J);   
    VecCreate(PETSC_COMM_WORLD,&x);
    VecSetSizes(x,PETSC_DECIDE,nodes.size());
    VecSetFromOptions(x);
    VecZeroEntries(x);
    //VecSet(x, 1.0);
    VecDuplicate(x,&r);
    VecDuplicate(x,&sol);
       
    array<int, 3> nodeIDs;
    array<array<double, 2>, 3> elemCoords;
    
    // initial guess
    /*for (vector<Element>::iterator it = getGrid().begin(); it != getGrid().end(); it++) {	
		nodeIDs[0] = it->getNode1().getID(); 
		nodeIDs[1] = it->getNode2().getID();
		nodeIDs[2] = it->getNode3().getID();
		elemCoords[0] = it->getNode1().getCoord();
		elemCoords[1] = it->getNode2().getCoord();
		elemCoords[2] = it->getNode3().getCoord();
		for (size_t i = 0; i < 3; i++) { // loop over nodes ('i'th equation/row for an element)			
			// ###### Jacobian matrix ######
			if ((elemCoords[i][0] != x_min) && (elemCoords[i][0] != x_max) && (elemCoords[i][1] != y_min) && (elemCoords[i][1] != y_max)) // left, right, bottom, top boundaries
				VecSetValues(x,1,&nodeIDs[i],&one,INSERT_VALUES);
		}
	}*/
				
	VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    //VecView(x,PETSC_VIEWER_STDOUT_SELF);		
}

// read mesh-data from file
void Mesh::readMesh() {
	ifstream mesh_read;
	mesh_read.open(meshFileName);
	// read node coordinates
	int n_nodes, n_dim; // number of nodes, number of dimensions	
	mesh_read >> n_nodes >> n_dim;
	nodes.resize(n_nodes);
	for (int i = 0; i < n_nodes; i++) {
		nodes[i].resize(n_dim);
	}
	for (int i = 0; i < n_nodes; i++) {
		for (int j = 0; j < n_dim; j++) {
			mesh_read >> nodes[i][j];
		}
	}
	
	// read connectivities
	int n_elems, n_nodes_per_elem;
	mesh_read >> n_elems >> n_nodes_per_elem; // number of elements, number of nodes per element
	elements.resize(n_elems);
	for (int i = 0; i < n_elems; i++) {
		elements[i].resize(n_nodes_per_elem);
	}
	for (int i = 0; i < n_elems; i++) {
		for (int j = 0; j < n_nodes_per_elem; j++) {
			mesh_read >> elements[i][j];
			elements[i][j] -= 1; // '0' - indexing
		}
	}
	mesh_read.close();

	for (int i = 0; i < n_nodes; i++) {
		if (x_max < nodes[i][0])
			x_max = nodes[i][0];
		if (x_min > nodes[i][0])
			x_min = nodes[i][0];
		if (y_max < nodes[i][1])
			y_max = nodes[i][1];
		if (y_min > nodes[i][1])
			y_min = nodes[i][1];
	}
}

// assemble global stiffness matrix
void Mesh::Assemble() {
	// loop through elements and add up contributions (except at Dirichlet boundary nodes)	
	MatZeroEntries(J);
	//PetscScalar one = 1.0;
	vector<int> nodeIDs;
	nodeIDs.resize(3);
	array<array<double, 2>, 3> elemCoords;
	int elemID;
	for (vector<Element>::iterator it = mesh.begin(); it != mesh.end(); it++) {
		MatAssemblyBegin(J,MAT_FLUSH_ASSEMBLY);
        MatAssemblyEnd(J,MAT_FLUSH_ASSEMBLY);
        elemID = it->getElemID();
		nodeIDs[0] = it->getNode1().getID(); 
		nodeIDs[1] = it->getNode2().getID();
		nodeIDs[2] = it->getNode3().getID();
		elemCoords[0] = it->getNode1().getCoord();
		elemCoords[1] = it->getNode2().getCoord();
		elemCoords[2] = it->getNode3().getCoord();              		
		const array<array<double, 3>, 3> elemStiffness = it->getElemStiffness();
		for (int i = 0; i < 3; i++) { 
			ierr = MatAssemblyBegin(J,MAT_FLUSH_ASSEMBLY);
			ierr = MatAssemblyEnd(J,MAT_FLUSH_ASSEMBLY);  
			if ((elemCoords[i][0] == x_min) || (elemCoords[i][0] == x_max) || (elemCoords[i][1] == y_min) || (elemCoords[i][1] == y_max)) { // left, right, bottom, top boundaries                        
				ierr = MatSetValues(J,1,&nodeIDs[i],1,&nodeIDs[i],&one,INSERT_VALUES);
			}
			else {
				ierr = MatAssemblyBegin(J,MAT_FLUSH_ASSEMBLY);
				ierr = MatAssemblyEnd(J,MAT_FLUSH_ASSEMBLY);                                
				for (int j = 0; j < 3; j++) {
					//globalStiffness[nodeIDs[i]][nodeIDs[j]] += elemStiffness[i][j];
					ierr = MatSetValues(J,1,&nodeIDs[i],1,&nodeIDs[j],&elemStiffness[i][j],ADD_VALUES);
					//cout << "Ag for " << nodeIDs[i] << ", " << nodeIDs[j] << ": " << globalStiffness[nodeIDs[i]][nodeIDs[j]] << endl;
				}				
			}
		}
	}
	//cout << "J assembled!" << endl;	
	MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);
    //MatView(J,PETSC_VIEWER_STDOUT_SELF);		
}

void Mesh::Solve() {
    PetscReal norm, tol = 1.e-14;      
    SNESCreate(PETSC_COMM_WORLD,&snes);
    SNESSetFunction(snes,this->r,FormFunction,this);    
    SNESSetJacobian(snes,this->J,this->J,FormJacobian,this);           
    SNESGetKSP(snes,&ksp);
    KSPGetPC(ksp,&pc);        
    PCSetType(pc,PCJACOBI);
    KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    double atol = 1e-12, rtol = 1e-12, stol = 1e-12;
    int maxit = 500, maxf = 1000;  
    SNESSetTolerances(snes, atol, rtol, stol, maxit, maxf);
    SNESSetFromOptions(snes);
    SNESSolve(snes,PETSC_NULL,x);
    copySol(x);
    SNESGetIterationNumber(snes,&its);
    //PetscPrintf(PETSC_COMM_WORLD,"number of SNES iterations = %D\n",its);
    //VecView(x,PETSC_VIEWER_STDOUT_SELF);   
}

// return reference to grid
vector<Element>& Mesh::getGrid() {
	return mesh;
}

Mat& Mesh::getJac() {
	//MatView(J,PETSC_VIEWER_STDOUT_SELF);     
    return J;
}
 
Vec& Mesh::getR() {
	return r;
}

const Vec& Mesh::getSol() {
	return x;
}
 
void Mesh::copySol(Vec& x_sol) {
	int n;
	//VecView(x_sol,PETSC_VIEWER_STDOUT_SELF);     
	//PetscScalar iVal;
	//VecGetSize(x_sol,&n);
	//Vec current_sol = getSol();
	/*for (int i = 0; i < n; i++) {
        VecGetValues(x_sol,1,&i,&iVal);
		VecSetValues(current_sol,1,&i,&iVal,ADD_VALUES);
	}*/
	//VecCopy(x_sol,current_sol);
	setSol(x_sol);	
	//VecView(current_sol,PETSC_VIEWER_STDOUT_SELF);     
} 
 
void Mesh::setSol(Vec& updated_sol) {
	//VecView(updated_sol,PETSC_VIEWER_STDOUT_SELF); 
	VecCopy(updated_sol,sol);
	//VecView(sol,PETSC_VIEWER_STDOUT_SELF); 
}

void Mesh::setR(Vec& r_iter) {
	VecCopy(r_iter,r);
//	VecView(r,PETSC_VIEWER_STDOUT_SELF); 
}

void Mesh::writeField(const string& solName) {
      ofstream outFile;           
      stringstream ss;
      ss << solName << ".plt";
      outFile.open(ss.str());
      outFile << "TITLE = \"SOLUTION_FIELD\"" << endl;
      outFile << "VARIABLES = \"x\" \"y\" \"u\"" << endl;
      outFile << "ZONE N=" << nodes.size() << ", E=" << elements.size() << ", F=FEPOINT, ET=TRIANGLE" << endl;
      x_sol.resize(nodes.size());
      for (int i=0; i<nodes.size(); i++) {
	      VecGetValues(x,1,&i,&x_sol(i));
              outFile << setprecision(6) << nodes[i][0] << "\t" << setprecision(6) << nodes[i][1] << "\t" << setprecision(6) << x_sol(i) << endl;
      }

      for (int i=0; i<elements.size(); i++) {
              outFile << elements[i][0] + 1 << "\t" << elements[i][1] + 1 << "\t" << elements[i][2] + 1 << endl;
      }
      outFile.close();        
}

// destructor
Mesh::~Mesh() { 
	MatDestroy(&J);
    VecDestroy(&x);  
    VecDestroy(&sol);
    SNESDestroy(&snes);
    PetscFinalize();
	cout << "Grid destroyed!" << endl; 
}


// residual calculation
PetscErrorCode FormFunction(SNES snes,Vec x,Vec f,void* meshPtr) {	
	PetscErrorCode ier;
    Mesh* mPtr = (Mesh*) meshPtr;
    /*cout << "r: " << endl;
    VecView(mPtr->r,PETSC_VIEWER_STDOUT_SELF);
    cout << "x: " << endl;
    VecView(mPtr->x,PETSC_VIEWER_STDOUT_SELF);
    cout << "J: " << endl;
    MatView(mPtr->getJac(),PETSC_VIEWER_STDOUT_SELF);*/    
    PetscScalar m_one = -1.0, zero = 0.0, a_val = 0.0;
    PetscScalar v, val1, val2, val3, v1, v2, v3, s1;
    double factor = 1/3.0;
	array<double, 3> gauss_pt_weights = { factor, factor, factor };
	Matrix<double, 3, 2> gauss_pts;
	gauss_pts << 0.5, 0.0, 0.0, 0.5, 0.5, 0.5;
	double evalPoint = 0.0, x1 = 0.0, x2 = 0.0, l = 0.0;
	array<double, 3> basis_values;
	Matrix<double, 3, 2> der_values;
	array<double, 2> coords;
	array<int, 3> nodeIDs;
	array<array<double, 2>, 3> elemCoords;
	mPtr->Assemble();
    VecZeroEntries(f);
	double exp_val;
	int elemID;
	double x_query, y_query;	
	//MatView(mPtr->J,PETSC_VIEWER_STDOUT_SELF); 
	array<array<double, 3>, 3> elemStiffness;
	for (vector<Element>::iterator it = mPtr->getGrid().begin(); it != mPtr->getGrid().end(); it++) {	
		elemStiffness = it->getElemStiffness();		
		nodeIDs[0] = it->getNode1().getID(); 
		nodeIDs[1] = it->getNode2().getID();
		nodeIDs[2] = it->getNode3().getID();
		elemCoords[0] = it->getNode1().getCoord();
		elemCoords[1] = it->getNode2().getCoord();
		elemCoords[2] = it->getNode3().getCoord();
	    VecGetValues(x,1,&nodeIDs[0],&val1);		        
		VecGetValues(x,1,&nodeIDs[1],&val2);
		VecGetValues(x,1,&nodeIDs[2],&val3);
		for (size_t i = 0; i < 3; i++) { // loop over nodes ('i'th equation/row for an element)			
			// ###### Jacobian matrix ######
			if ((elemCoords[i][0] != mPtr->x_min) && (elemCoords[i][0] != mPtr->x_max) && (elemCoords[i][1] != mPtr->y_min) && (elemCoords[i][1] != mPtr->y_max)) { // left, right, bottom, top boundaries 
				for (size_t j = 0; j < 3; j++) { // 'j'th column for 'i'th row
					v = 0.0;
					for (size_t k = 0; k < gauss_pts.rows(); k++) {
						basis_values = it->getBasis().calcBasis(gauss_pts(k, 0), gauss_pts(k, 1));
						exp_val = val1*basis_values[0] + val2*basis_values[1] + val3*basis_values[2]; // hard-coded for linear triangles				
						v += 0.5 * gauss_pt_weights[k] * basis_values[i] * basis_values[j] * pow(PI, 2.0) * exp(exp_val) * it->getBasis().getDetJ();					
					}
					MatSetValues(mPtr->J,1,&nodeIDs[i],1,&nodeIDs[j],&v,ADD_VALUES);	
				}	    
		    //MatAssemblyBegin(mPtr->getJac(),MAT_FLUSH_ASSEMBLY);
            //MatAssemblyEnd(mPtr->getJac(),MAT_FLUSH_ASSEMBLY);		
	        }
	        
			// ###### residual vector ######			
			v = 0.0; v1 = 0.0; v2 = 0.0; v3 = 0.0; s1 = 0.0;		
			if ((elemCoords[i][0] != mPtr->x_min) && (elemCoords[i][0] != mPtr->x_max) && (elemCoords[i][1] != mPtr->y_min) && (elemCoords[i][1] != mPtr->y_max)) { // left, right, bottom, top boundaries 		
				der_values = it->getBasis().calcBasisDer(); // should take in Gauss point for higher-order basis functions, but here it is outside loop because of linear basis - WARNING!!!!
				v1 = (pow(der_values(0,0), 2.0) + pow(der_values(0,1), 2.0)) * val1 + (der_values(0,0) * der_values(1,0) + der_values(0,1) * der_values(1,1)) * val2 + (der_values(0,0) * der_values(2,0) + der_values(0,1) * der_values(2,1)) * val3;
				v2 = (pow(der_values(1,0), 2.0) + pow(der_values(1,1), 2.0)) * val2 + (der_values(0,0) * der_values(1,0) + der_values(0,1) * der_values(1,1)) * val1 + (der_values(1,0) * der_values(2,0) + der_values(1,1) * der_values(2,1)) * val3;
				v3 = (pow(der_values(2,0), 2.0) + pow(der_values(2,1), 2.0)) * val3 + (der_values(0,0) * der_values(2,0) + der_values(0,1) * der_values(2,1)) * val1 + (der_values(1,0) * der_values(2,0) + der_values(1,1) * der_values(2,1)) * val2;
	            for (size_t k = 0; k < gauss_pts.rows(); k++) {
		            basis_values = it->getBasis().calcBasis(gauss_pts(k, 0), gauss_pts(k, 1));
					exp_val = val1*basis_values[0] + val2*basis_values[1] + val3*basis_values[2];
					v += 0.5 * gauss_pt_weights[k] * basis_values[i] * pow(PI, 2.0) * exp(exp_val) * it->getBasis().getDetJ();
					s1 += 0.5 * gauss_pt_weights[k] * it->getBasis().getDetJ();
	            }	
	            if (i == 0) {
					v1 *= s1;           
					v += v1;
				}
				else if (i == 1) {
					v2 *= s1;           
					v += v2;
				}
				
				else if (i == 2) {
					v3 *= s1;           
					v += v3;
				}
				ier = VecSetValues(f,1,&nodeIDs[i],&v,ADD_VALUES); CHKERRQ(ier);
				//VecView(f,PETSC_VIEWER_STDOUT_SELF); 		    
	        }
	    }
    }    
    //VecCopy(f,mPtr->r);
    //cout << "iter: " << iter << endl;
    //VecView(f,PETSC_VIEWER_STDOUT_SELF); 
    //mPtr->copySol(x);
    //mPtr->setR(f);
    //VecAssemblyBegin(mPtr->getR());
    //VecAssemblyEnd(mPtr->getR());
    MatAssemblyBegin(mPtr->J,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mPtr->J,MAT_FINAL_ASSEMBLY);	     
    //MatView(mPtr->J,PETSC_VIEWER_STDOUT_SELF);
  //  VecView(mPtr->getR(),PETSC_VIEWER_STDOUT_SELF);
    //iter++;
}

// Jacobian 
PetscErrorCode FormJacobian(SNES snes,Vec x,Mat jac,Mat B,void* meshPtr){	
	return 0;   
}
