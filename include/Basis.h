#pragma once

#include"Node.h"
#include <array>
#include<Eigen/Dense>

class Basis {
public:
	Basis(const Node&, const Node&, const Node&); // default constructor
	std::array<double, 3> calcBasis(const double&, const double&); // returns 3 x 1 vector of basis function values at query point (for triangle elements)
	Eigen::Matrix<double, 3, 2> calcBasisDer(); // returns 3 x 2 dNdX matrix of basis function derivative values at query point (for triangle elements)
	void calcDEDX(); // return dEdN/dXdY
	void calcDetJ(); // calculate determinant of Jacobian
	const double getDetJ(); // return determinant of Jacobian
        ~Basis(); // destructor

private:
	Eigen::Matrix<double, 3, 2> dNdE = Eigen::ArrayXXd::Zero(3,2);
	Eigen::Matrix<double, 2, 2> dEdX = Eigen::ArrayXXd::Zero(2,2);
	std::array<std::array<double, 2>, 3> coords;
	double detJ;
};
