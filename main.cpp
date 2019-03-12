#include <iostream>
#include <fstream>
#include "ele.h"
#include "reconstruct.h"
#include "quadrature.h"
#include "coordtransform.h"
Polynomial pp;
int main(int argc, char ** argv) {
  double left, right,erro;
  int n_ele, n_patch, order;

  if(argc < 6) {
    std::cerr << "not enough parameters" << std::endl;
    return 1;
  }

  left = atof(argv[1]);  right = atof(argv[2]);
  n_ele = atoi(argv[3]); n_patch = atoi(argv[4]);
  order = atoi(argv[5]);

  std::vector<Element> ele_cache;

  pp.set_n_array(order);

  reconstruct(ele_cache, left, right, n_ele, n_patch, order);
  
  int quad_acc = 8;

  find_quadrature_info(ele_cache, left, right, n_ele, quad_acc);
  
  erro = L2_erro(ele_cache); 

  std::ofstream output2("l2erro.txt");
  output2<<erro<<std::endl;
  output2.close();

}
