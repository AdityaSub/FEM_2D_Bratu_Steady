#pragma once
#include<vector>
#include<Eigen/Dense>
#include "Element.h"
#include<petscsnes.h>

class Mesh {
public:
	Mesh(std::string&); // constructor
	void readMesh(); // read mesh-data from file
	std::vector<Element>& getGrid(); // return reference to grid
	void Assemble(); // assemble linear system
	void Solve(); // solve linear system for unknowns
	Mat& getJac(); // return Jacobian-matrix
	Vec& getR();
	const Vec& getSol(); // return solution-vector
	void copySol(Vec&);
	void setSol(Vec&);
	void setR(Vec&);
	void writeField(const std::string&);
	Eigen::VectorXd& getNodeList();
	void computeL2Error(); // calculate L2 error
	~Mesh(); // destructor
	Vec x, r;
    Mat J;
    double x_min = 0.0, y_min = 0.0, x_max = 0.0, y_max = 0.0; // mesh bounds

private:
	std::string meshFileName; // input mesh-filename
	std::vector <std::vector<double>> nodes; // node-list for grid
	std::vector <std::vector<int>> elements; // connectivities
	std::vector<Element> mesh; // grid-vector containing list of elements
	std::vector<double> analyticalSolution; // analytical solution vector	
	PetscScalar one = 1.0, zero = 0.0;	
	Vec sol;
    SNES snes;
    KSP ksp;
    PC pc;
    PetscErrorCode ierr;
    PetscInt its;
    Eigen::VectorXd x_sol; // solution-vector
};
