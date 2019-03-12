#include <vector>
#include "Eigen/Dense"
#include <math.h>
#include <iostream>
#include <algorithm>
Eigen::MatrixXd get_weight(
    std::vector<double> x,
    std::vector<int> free,
    double x0,
    int m)
{
    Eigen::MatrixXd WW = Eigen::MatrixXd::Identity(free.size(),free.size());
    Eigen::MatrixXd CC = Eigen::MatrixXd::Zero(m+1,m+1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(m+1);
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(m+1,free.size());
    std::vector<int> free_cache = free;
    std::vector<int> freevar,notfreevar;
    int cindex;
    std::vector<double> a(m+1);
    switch(m)
    {
        case 2: 
        a[0] = 0; a[1] = 10; a[2] = 1;
        if (free.size()<9) return WW;
        break;
    }
   // for(auto& sb:free)
    //    std::cout<<sb<<std::endl;
    for(int i = 0; i<free.size(); i++)
        {
            if(x[free[i]]==x0)
             {
              freevar.push_back(i);
              cindex = i;
              break;
             }
        }
    double h = fabs(x[free[1]]-x[free[0]]);
    //std::cout<< a[0] <<" "<<a[1]<<" "<<a[2]<< std::endl;
    std::vector<double> c;
    for(int i = 0; i<free.size(); i++)
        c.push_back(free[i]-free[cindex]);
    //for(auto& sb:c)
    //    std::cout<<sb<<std::endl;
    for(int j = 0; j<m+1; j++)
    {
        for(int i = 0; i<free.size(); i++)
        {
            for (int k = 0; k<m+1; k++)
                temp(j,i) = temp(j,i)+a[k]*pow(c[i],k);
            temp(j,i) = temp(j,i)-pow(c[i],m+1);
            temp(j,i) = temp(j,i)*pow(h,m+1+j)*pow(c[i],j);
            temp(j,i) = temp(j,i)*pow(10,7);
        }
    }
   //std::cout <<"This is temp"<<std::endl;
   //std::cout<<temp<<std::endl; 
   /* bool neg = false; bool pos = false;
    
    for(int i = 0; i<free.size(); i++)
    {
        if (i == cindex) continue;
        if ((neg == false)&&(notfreevar.size()<m+1)&&(temp(0,i)<0))
        {
            notfreevar.push_back(i);
            neg = true;
        }
        else
        if ((pos == false)&&(notfreevar.size()<m+1)&&(temp(0,i)>0))
        {
            notfreevar.push_back(i);
            pos = true;
        }
    }
    //std::cout<<"hahaha"<<std::endl;
    for(int i = 0; i<free.size(); i++)
    {
        if((find(notfreevar.begin(),notfreevar.end(),i) == notfreevar.end())&&(find(freevar.begin(),freevar.end(),i) == freevar.end()))
        {
            if(notfreevar.size()<m+1) notfreevar.push_back(i); else freevar.push_back(i);    
        }
    }
   */
   for(int i = free.size()-1; i>=0; i--)
    {
        if(x[free[i]]==x0) continue;
        if(notfreevar.size()<m+1) notfreevar.push_back(i);
            else freevar.push_back(i);
    }
    sort(notfreevar.begin(),notfreevar.end());
     for(int i = 1; i<freevar.size(); i++)
         WW(freevar[i],freevar[i]) = 1;
         WW(freevar[3],freevar[3]) = 10.5;
   // for(auto& sb:freevar)
   //     std::cout<<sb<<" ";
   // std::cout<<std::endl;
    for(int j = 0; j<m+1; j++)
    {
        double sum = 0;
        for(int i = 1; i<freevar.size(); i++)
        {
            sum = sum-WW(freevar[i],freevar[i])*temp(j,freevar[i]);
        }
        b(j) = sum;
        //std::cout<<b(j)<<" ";
        for(int i = 0; i<notfreevar.size(); i++)
        {
            CC(j,i) = temp(j,notfreevar[i]);
        }
    }
   // std::cout<<CC<<std::endl; 
    //std::cout<<std::endl;
   // std::cout<<b<<std::endl;
    Eigen::VectorXd s = CC.fullPivLu().solve(b);
    WW(freevar[0],freevar[0]) = 1000; //中心的权重大一些 
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
