#include<iostream>
#include "Node.h"

using namespace std;

// constructor w/ initializer list
Node::Node(int ID, double loc_x, double loc_y) :nodeID(ID), coord({ loc_x, loc_y }) {}

// copy constructor
Node::Node(const Node& N) {
	nodeID = N.nodeID;
	coord[0] = N.coord[0];
	coord[1] = N.coord[1];
}

// return location of current node
const array<double, 2>& Node::getCoord() const {
	return coord;
}

// return ID of current node
const int& Node::getID() const {
	return nodeID;
}

// destroy current node object
Node::~Node(){}