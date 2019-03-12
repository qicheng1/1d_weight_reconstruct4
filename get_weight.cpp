#include <vector>
#include "Eigen/Dense"
#include <math.h>
#include <iostream>
Eigen::MatrixXd get_weight(
    std::vector<double> x,
    std::vector<int> free,
    double x0,
    double m)
{
    Eigen::MatrixXd WW = Eigen::MatrixXd::Zero(free.size(),free.size());
    Eigen::MatrixXd CC = Eigen::MatrixXd::Zero(m+1,m+1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(m+1);
    std::vector<int> freevar,notfreevar;
    for(int i = 0; i<free.size(); i++)
        {
            if(x[free[i]]==x0)
             {
              freevar.push_back(i);
              break;
             }
        }
    int cindex = i;
    double h = fabs(x[free[1]]-x[free[0]]);
    std::vector<double> a(m+1);
    if (m%2==0) a[0]=0; else a[0]=pow(10,m-1);
    for(int i = 1; i<m+1; i++)
    {
        a[i]=pow(10,m-i-1);
    }
    std::vector<double> c;
    for(int i = 0; i<free.size(); i++)
    {
        if(x[free[i]]==x0)
         {
             c.push_back(free[i]-free[cindex]);
             continue;
         }
         c.push_back(free[i]-free[cindex]);
        if(notfreevar.size()<m+1) notfreevar.push_back(i);
        else freevar.push_back(i);
    }

    for(int j = 0; j<m+1; j++)
    {
        double sum = 0;
        for(int i = 1; i<freevar.size(); i++)
        {
            double temp;
            for (int k = 0; k<m+1; k++)
                temp = temp+a[k]*pow(c[freevar[i]],k);
            temp = temp-pow(c[freevar[i]],m+1);
            temp = temp*pow(h,m+1+j)*pow(c[freevar[i]],j);
            sum = sum-1*temp;
        }
        b(j) = sum;
        for(int i = 0; i<notfreevar.size(); i++)
        {
            double temp1;
            for (int k = 0; k<m+1; k++)
                temp1 = temp1+a[k]*pow(c[notfreevar[i]],k);
            temp1 = temp1-pow(c[notfreevar[i]],m+1);
            temp1 = temp1*pow(h,m+1+j)*pow(c[notfreevar[i]],j);
            CC(j,i) = temp1;
        }
    } 
    Eigen::VectorXd s = CC.fullPivLu().solve(b);
    WW(freevar[0],freevar[0]) = 1000; //中心的权重大一些
    for(int i = 1; i<freevar.size(); i++)
        WW(freevar[i],freevar[i]) = 1; 
    for(int i = 0; i<notfreevar.size(); i++)
        WW(notfreevar[i],notfreevar[i]) = s[i];
    for(int i = 0; i<free.size(); i++)
       std::cout<<WW(i,i)<<std::endl;
    std::cout<<std::endl;
    //Eigen::MatrixXd WW = Eigen::MatrixXd::Identity(free.size(),free.size());
/*    std::vector<double> distance;
     double s = 0;
    for(auto &i : free)
    {
	distance.push_back(fabs(x[i]-x0));
	if(fabs(x[i]-x0) != 0) s+=1.0/fabs(x[i]-x0);
    }
    for(auto k = 0; k<distance.size(); k++)
    {
	if(distance[k]==0) distance[k] = 1000; //单元中心的权近似无穷大
	else
	    distance[k] = 1.0/(distance[k]*s); //其他点的权重为距离倒数（归1化）
    }

    //std::cout<<"s = "<<s<<std::endl;
    
    for(auto j = 0; j<free.size(); j++)
	    WW(j,j) = distance[j];
*/
    
    return WW;	
}
