#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;


int main() {

	Matrix2d a;
	a << 1,2,3,4;

	std::cout << a << std::endl;
}