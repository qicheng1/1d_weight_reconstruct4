#ifndef __QUADRATURE_H_
#define __QUADRATURE_H_

#include <vector>
#include <cassert>

class QuadratureInfo{
public:
  typedef double point_t;
public:
  QuadratureInfo() = default; //Default contructor//
  QuadratureInfo(const int & a,const std::vector< point_t >& vec1,const std::vector<double>& vec2):acc(a),point(vec1),wei(vec2) {}; //Initial contructor//
  QuadratureInfo(const QuadratureInfo& quad):acc(quad.algebricAccuracy()),point(quad.quadraturePoint()),wei(quad.weight()) {}; //Copy contructor//
  ~QuadratureInfo() = default; //Destructor//
public:
  QuadratureInfo& operator=(const QuadratureInfo&); //Copy operator//
  int n_quadraturePoint() const; // Number of quadrature points//
  int algebricAccuracy() const; //Algebraic accuracy//
  int &algebricAccuracy(); //Algebraic accuracy//
  const std::vector< point_t >& quadraturePoint() const; // Quadrature point array//
  std::vector< point_t >& quadraturePoint(); // Quadrature point array//
  const point_t& quadraturePoint(int) const; // Quadrature point//
  point_t& quadraturePoint(int); // Quadrature point//
  const std::vector<double>& weight() const; // Quadrature weight array //
  std::vector<double>& weight(); // Quadrature weight array //
  const double& weight(int) const; //Quadrature weight//
  double& weight(int); //Quadrature weight//

private:
  int acc;
  std::vector<point_t> point;
  std::vector<double> wei;
};


class TemplateQuadratureInfo{
public:
  typedef QuadratureInfo quad_t;
  typedef double point_t;

  quad_t& findQuadratureInfo(int i){
    assert(i <= acc_tbl[acc_tbl.size() - 1]);
    if(i == acc_tbl[0]) return quad[0];

    for(int j = 0; j < acc_tbl.size() - 1; ++j)
      if(i <= acc_tbl[j + 1] && i > acc_tbl[j])
        return quad[j + 1];
  }

  const quad_t& findQuadratureInfo(int i) const {
    assert(i <= acc_tbl[acc_tbl.size() - 1]);
    if(i == acc_tbl[0]) return quad[0];

    for(int j = 0; j < acc_tbl.size() - 1; ++j)
      if(i <= acc_tbl[j + 1] && i > acc_tbl[j])
        return quad[j + 1];
  }

protected:
  std::vector<int> acc_tbl; // Algebraic accuracy table. //
  std::vector<quad_t> quad; // QuadratureInfo table//
};

class LineQuadratureInfo : public TemplateQuadratureInfo{
  typedef TemplateQuadratureInfo Base;
  using Base::acc_tbl;
  using Base::quad;
  typedef double point_t;
public:
  LineQuadratureInfo() {
  acc_tbl={1,3,7,11};
  quad={
      QuadratureInfo(1,
        {point_t(0.0)},
        {1.0}),
      QuadratureInfo(3,
        {point_t(-0.57735026918963),point_t(0.57735026918963)},
        {0.5,0.5}),
      QuadratureInfo(7,
        {point_t(0.211324865),point_t(0.788675135),point_t(-0.211324865),point_t(-0.788675135)},
        {0.25,0.25,0.25,0.25}),
      QuadratureInfo(11,
        {point_t(0.112701665),point_t(0.500000000),point_t(0.887298335),point_t(-0.112701665),point_t(-0.500000000),point_t(-0.887298335)},
        {0.138888889,0.222222222,0.138888889,0.138888889,0.222222222,0.138888889})
      };
  }
};


#endif
