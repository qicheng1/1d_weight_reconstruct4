#include <vector>
#include "ele.h"
#include "math.h"
#include <iostream>
double g(double x);
double L2_erro(std::vector<Element>& ele)
{
	int i,j,k,l;
	double temp,ans;
	ans = 0;
	for(l = 0; l<ele.size(); l++)
	{
		std::vector<double> nodes(ele[l].get_interval_cache().quad.quadraturePoint());
		std::vector<double> wei(ele[l].get_interval_cache().quad.weight());
		double jacobi = ele[l].get_interval_cache().jacobian[0];
		std::vector<int> free(ele[l].get_dof());
		for(i = 0; i<nodes.size(); i++)
		{
			temp = 0;
			std::vector<double> basis_value(ele[l].basis_function_value(nodes[i]));
			for(k=0; k<ele[l].n_dof(); k++)
			  temp = temp+basis_value[k]*g(ele[free[k]].get_center());
			temp = temp-g(nodes[i]);
           // std::cout<< temp/pow(nodes[i]-ele[l].get_center(),3) << " ";
			ans = ans+temp*temp*wei[i]*jacobi;
		}
        //std::cout<<std::endl;
	}
	return sqrt(ans);
}
