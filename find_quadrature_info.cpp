#include "ele.h"
#include <vector>
#include "quadrature.h"
#include "coordtransform.h"
void find_quadrature_info(std::vector<Element>& ele, 
    double a, 
    double b, 
    int n, 
    int acc)
{
    int i,j,k,l;
    std::vector<double> xx;
    std::vector<double> lv = {-1,1};
    LineQuadratureInfo quadifo;
    CoordTransform coordtrans;
    QuadratureInfo quad = quadifo.findQuadratureInfo(acc);
    for(i = 0; i<=n; i++)
    {
		xx.push_back(a+(b-a)*i/n);
    }
    for(l = 0; l<ele.size(); l++)
    {
    	std::vector<double> nodes(quad.quadraturePoint());
    	std::vector<double> gv = {xx[l],xx[l+1]};
        for(i = 0; i<nodes.size(); i++)
        {
        	nodes[i] = coordtrans.local_to_global(nodes[i],lv,gv);
        }
        std::vector<double> jaco = {coordtrans.local_to_global_jacobian(0,lv,gv)};
        ele[l].set_cache(xx[l],xx[l+1],QuadratureInfo(quad.algebricAccuracy(),nodes,quad.weight()),jaco);
	}
}
