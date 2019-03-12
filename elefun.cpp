#include <vector>
#include "ele.h"
#include "math.h"
#include "quadrature.h"
void Element::set_cache(double a, double b, QuadratureInfo qu, std::vector<double> ja)
{
	cache.left = a;
        cache.right = b;
        cache.quad = qu;
        cache.jacobian = ja;
}
std::vector<double> Element::basis_function_value(const double x) const
{
	std::vector<double> vec(dof.size());
	for(auto i = 0; i<dof.size(); i++)
		{
			vec[i] = poly->value(x,basis[i],center);
		}
	return vec;
}
std::vector<std::vector<double>> Element::basis_function_gradient(const double x) const
{
	std::vector<std::vector<double>> mtx(dof.size(),std::vector<double>(poly->get_order()));
	for(auto i = 0; i<dof.size(); i++)
		{
			mtx[i] = poly->gradient(x,basis[i],center);
		}
	return mtx;
}
double Polynomial::value(const double x, const std::vector<double>& coe,const double x0 = 0) const
{
    double ans = 0;
    for(auto k = 0; k<coe.size(); k++)
		ans = ans+pow(x-x0,k)*coe[k];
	return ans;
}
std::vector<double> Polynomial::gradient(const double x, const std::vector<double>& coe,const double x0 = 0) const
{
	std::vector<double> ans(order);
	for(auto k = 1; k<coe.size(); k++)
		ans[0] = ans[0]+k*pow(x-x0,k-1)*coe[k];
	return ans;
}

