#pragma once
#include<array>

class Node {
public:
	Node(const int = 0, const double = 0.0, const double = 0.0); // constructor
	Node(const Node &); // copy constructor
	const std::array<double, 2>& getCoord() const; // return coordinates of current node 
	const int& getID() const; // return node ID of current node
	~Node(); // destructor

private:
	int nodeID; // node ID of current node
	std::array<double, 2> coord; // location of current node - (x, y)
};