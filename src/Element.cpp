#include<iostream>
#include<iomanip>
#include<vector>
#include<math.h>
#include "Element.h"
#include "Basis.h"
#include<Eigen/Dense>

using namespace std;
using namespace Eigen;

// constructor 
Element::Element(const Node& N1, const Node& N2, const Node& N3, const int& ID)
	:node1(N1), node2(N2), node3(N3), elemID(ID), basis(Basis(node1, node2, node3)) {
        for (size_t i = 0; i < elemStiffness.size(); i++) {
                for (size_t j = 0; j < elemStiffness[0].size(); j++)
                        elemStiffness[i][j] = 0.0;
        }
        //cout << "Initialized element: " << elemID << endl;
}

// return 'node1'
const Node& Element::getNode1() const{
	return node1;
}

// return 'node2'
const Node& Element::getNode2() const{
	return node2;
}

// return 'node3'
const Node& Element::getNode3() const {
	return node3;
}

// return element ID
const int& Element::getElemID() const {
	return elemID;
}

// return basis-object for current element
Basis& Element::getBasis() {
        return basis;
} 

// calculate element stiffness matrix
void Element::calcElementStiffness() {
	double factor = 1/3.0;
	array<double, 3> gauss_pt_weights = { factor, factor, factor };
	Matrix<double, 3, 2> gauss_pts;
	gauss_pts << 0.5, 0.0, 0.0, 0.5, 0.5, 0.5;
	Matrix<double, 3, 2> der_values = basis.calcBasisDer(); // constant for linear elements
	for (int i = 0; i < gauss_pts.rows(); i++) {
		elemStiffness[0][0] += (0.5) * gauss_pt_weights[i] * (pow(der_values(0, 0), 2.0) + pow(der_values(0, 1), 2.0)) * basis.getDetJ();
		elemStiffness[0][1] += (0.5) * gauss_pt_weights[i] * (der_values(0, 0) * der_values(1, 0) + der_values(0, 1) * der_values(1, 1)) * basis.getDetJ();
		elemStiffness[0][2] += (0.5) * gauss_pt_weights[i] * (der_values(0, 0) * der_values(2, 0) + der_values(0, 1) * der_values(2, 1)) * basis.getDetJ();
		elemStiffness[1][1] += (0.5) * gauss_pt_weights[i] * (pow(der_values(1, 0), 2.0) + pow(der_values(1, 1), 2.0)) * basis.getDetJ();
		elemStiffness[1][2] += (0.5) * gauss_pt_weights[i] * (der_values(1, 0) * der_values(2, 0) + der_values(1, 1) * der_values(2, 1)) * basis.getDetJ();
		elemStiffness[2][2] += (0.5) * gauss_pt_weights[i] * (pow(der_values(2, 0), 2.0) + pow(der_values(2, 1), 2.0)) * basis.getDetJ();		
	}
	elemStiffness[1][0] = elemStiffness[0][1]; // symmetric K
	elemStiffness[2][0] = elemStiffness[0][2];
	elemStiffness[2][1] = elemStiffness[1][2];      
    
        //cout << "K_elem computed for elem: " << elemID << endl;
}

const array<array<double, 3>, 3> Element::getElemStiffness() const {
	return elemStiffness;
}

void Element::printElemStiffness() const {
        for (size_t i = 0; i < elemStiffness.size(); i++) {
                for (size_t j = 0; j < elemStiffness[0].size(); j++)
                        cout << setprecision(6) << setw(20) << elemStiffness[i][j];
                cout << endl;
        }
        cout << endl;
}

// destructor
Element::~Element() {}
