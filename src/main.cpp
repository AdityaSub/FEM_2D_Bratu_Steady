// main.cpp
// Solve a 2D heat-conduction problem with Dirichlet boundary conditions (u = 1 on top face, u = 0 on sides and bottom face)

#include<iostream>
#include<vector>
#include<string>
#include "Mesh.h"

using namespace std;

static char help[] = "Solves 2D Poisson's equation with KSP.\n\n";

int main(int argc, char* argv[])
{
	PetscInitialize(&argc,&argv,(char*)0,help);
    
	string fileName = argv[1];
	Mesh m(fileName); 
	//m.Assemble(); 
	//m.writeField("initial"); 
	m.Solve(); 
	m.writeField("solution");
}

