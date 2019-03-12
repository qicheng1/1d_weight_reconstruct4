#ifndef __ELE_H_
#define __ELE_H_

#include <vector>
#include "quadrature.h"

struct IntervalCache {
  double left, right;
  QuadratureInfo quad;
  std::vector<double> jacobian;
};


class Polynomial {
public:
  Polynomial() : order(1), n_array(2) {};
  Polynomial(int n): order(1), n_array(n){};

  void set_order(int d) { order = d;};
  void set_n_array(int n) { n_array = n;};
  int get_order() const {return order;};
  int get_n_array() const {return n_array;};

  double value(const double x, const std::vector<double>& coe, const double x0) const;
  //coe={a, b, c, ...} 代表多项式a + bx + cx^2 + ... 返回这个多项式在x处的值

  std::vector<double> gradient(const double x, const std::vector<double>& coe, const double x0) const;

  //void basis_array(const double x, std::vector<double>& array) const;
  //生成1 x x^2 x^3 ... 

private:
  int order;
  int n_array;
};

class Element {
public:
  Element (std::vector<int> vec, 
      double x,
      Polynomial & p, 
      std::vector<std::vector<double>> mtx):
    dof(vec),
    center(x),
    poly(&p),
    basis(mtx) {};
  // mtx存一个(A^TA)^-1A^T的转置
  Element() = default;
  // ... 构造函数
  std::vector<double> basis_function_value(const double x) const; 
  // 在x处基函数的值

  std::vector<std::vector<double>> basis_function_gradient(const double x) const; 
  //在x处插值多项式的梯度
  
  int n_dof() const { return dof.size();}

  std::vector<int>& get_dof() { return dof;}
  const std::vector<int>& get_dof() const { return dof;}
  
  //int get_center() const { return center;}
  double get_center() const { return center;}
  
  IntervalCache& get_interval_cache() { return cache;}
  const IntervalCache& get_interval_cache() const { return cache;}

  void set_cache(double a, double b, QuadratureInfo qu, std::vector<double> ja);
  // ... 

private:
  Polynomial * poly;
  int index; //第几个单元
  double center; //中心坐标
  std::vector<int> dof; // 这个单元上有几个自由度 也就是patch中的单元编号
  std::vector<std::vector<double>> basis; 
  //基函数 实际上就是(A^TA)^-1A^T basis[i]应该是第i个自由度对应的基函数
  //比如 a + bx + cx^2 那么basis[i]应该是{a, b, c}
  
  IntervalCache cache;

};

#endif
