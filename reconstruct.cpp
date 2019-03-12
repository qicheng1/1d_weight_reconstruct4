#include "ele.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "reconstruct.h"
#include "math.h"
#include "Eigen/Dense"
extern Polynomial pp;
void reconstruct(std::vector<Element>& ele, double a, double b, int n, int n_p, int m)
{
	std::ofstream output("basisvector.txt");
	std::vector<double> x;
	int i,j,k,l;
	for(i = 0; i<n; i++)
	{
		x.push_back(a+(b-a)/(2*n)+(b-a)*i/n);
	}
	//double erro_sum = 0; //总偏差平方和
    for(l = 0; l<n; l++) //处理第l个模板
    {  	
        std::vector<std::vector<double>> basis;
        std::vector<int> free;
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(m+1,m+1);
         Eigen::MatrixXd L_T = Eigen::MatrixXd::Zero(m+1,m+1);
         Eigen::MatrixXd U_T = Eigen::MatrixXd::Zero(m+1,m+1);
         Eigen::MatrixXd L_Ti = Eigen::MatrixXd::Zero(m+1,m+1);
         Eigen::MatrixXd U_Ti = Eigen::MatrixXd::Zero(m+1,m+1);
        //求第l个模板自由度
    	for(k = l-n_p; k<=l+n_p; k++)
    	{
    		if (k<0 || k>=n) continue;
    		free.push_back(k);
    	}
    	Eigen::MatrixXd T1 = Eigen::MatrixXd::Zero(m+1,free.size());
    	//计算系数矩阵的转置，保存在T1中
    	for(j = 0; j<m+1; j++)
    	    for(i = 0; i<free.size(); i++)
    	        T1(j,i) = pow(x[free[i]]-x[l],j);
	//std::cout << "The T1 is :\n" << T1 <<std::endl;
	Eigen::MatrixXd W = get_weight(x,free,x[l],m);
        //std::cout << "The W is :\n" << W <<std::endl; 
    	//计算(T1*W*'T1)
	   T = T1*W;
	   T = T*T1.transpose();
	   //std::cout << "The T is :\n" << T <<std::endl;
	   //lu_decompose(T,L_T,U_T,T1);
	   //求L_T,U_T的逆矩阵
	   //lu_inverse(L_T,U_T,L_Ti,U_Ti);
	   //计算基函数，即(A^TA)^-1A^T的转置
	   //std::cout<< "The LUparts of T \n" << T.fullPivLu().matrixLU() <<std::endl;
           Eigen::MatrixXd Ti = T.fullPivLu().inverse();
	   //std::cout << "The T*T-1 is \n" << T*Ti <<std::endl;
	   basis = get_basis(Ti,T1,W);
         // for (auto &e:basis)
         //     e.insert(e.begin(),0); 
	   output<<"The basis of element "<<l<<" :"<<std::endl;
	  for(i = 0; i<basis.size(); i++)
	  {
	      for(j = 0; j<basis[0].size(); j++)
		  output<<basis[i][j]<<" ";
	      output<<std::endl;
	  }
	  output<<std::endl;
	  ele.emplace_back(free,x[l],pp,basis);
    }
	
	/*std::vector<double> ba_v;
	erro_sum = 0; double t1 = 0;
	for(i = 0; i<n; i++)
	    {
	    	ba_v = elem[i].basis_function_value(x[i]);
	    	t1 = 0;
	    	for(k = 0; k<elem[i].n_dof(); k++)
	    	    t1 = t1+ba_v[k]*g(x[elem[i].get_dof()[k]]);
                \\std::cout<<t1<<std::endl;
	    	erro_sum = erro_sum+(t1-g(x[i]))*(t1-g(x[i]));
	    }
	output<<"erro_sum = "<<erro_sum<<std::endl;*/
	output.close();
}

