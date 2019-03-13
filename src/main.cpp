// main.cpp
// Solve a 2D Bratu equation (del^2 u = pi^2*exp(u)) with Dirichlet boundary conditions (u = 0 on all sides)

#include<iostream>
#include<vector>
#include<string>
#include "Mesh.h"

using namespace std;

static char help[] = "Solves 2D Bratu equation with SNES.\n\n";

int main(int argc, char* argv[])
{
	PetscInitialize(&argc,&argv,(char*)0,help);
    
	string fileName = argv[1]; // mesh file-name
	Mesh m(fileName);  // read mesh-data: nodes, connectivities
	//m.Assemble();
	//m.writeField("initial"); // write initial field to file
	m.Solve();  // solve non-linear functional
	m.writeField("solution"); // write solution data to file
}

