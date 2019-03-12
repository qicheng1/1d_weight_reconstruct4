#include "quadrature.h"
#include <vector>
QuadratureInfo& QuadratureInfo::operator=(const QuadratureInfo& quad) //Copy operator//
{
	acc = quad.algebricAccuracy();
	point = quad.quadraturePoint();
	wei = quad.weight();
	return *this;
}
int QuadratureInfo::n_quadraturePoint() const // Number of quadrature points//
{
    return point.size();
}
int QuadratureInfo::algebricAccuracy() const //Algebraic accuracy//
{
	return acc;
}
int &QuadratureInfo::algebricAccuracy() //Algebraic accuracy//
{
    return acc;
}
const std::vector< QuadratureInfo::point_t >& QuadratureInfo::quadraturePoint() const // Quadrature point array//
{
	return point;
}
std::vector< QuadratureInfo::point_t >& QuadratureInfo::quadraturePoint() // Quadrature point array//
{
	return point;
}
const QuadratureInfo::point_t& QuadratureInfo::quadraturePoint(int i) const // Quadrature point//
{
	return point[i];
}
QuadratureInfo::point_t& QuadratureInfo::quadraturePoint(int i) // Quadrature point//
{
	return point[i];
}
const std::vector<double>& QuadratureInfo::weight() const // Quadrature weight array //
{
    return wei;
}
std::vector<double>& QuadratureInfo::weight() // Quadrature weight array //
{
	return wei;
}
const double& QuadratureInfo::weight(int i) const //Quadrature weight//
{
    return wei[i];
}
double& QuadratureInfo::weight(int i) //Quadrature weight//
{
    return wei[i];
}
