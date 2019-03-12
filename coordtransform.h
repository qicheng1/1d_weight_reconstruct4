#ifndef __COORDTRANSFORM_H_
#define __COORDTRANSFORM_H_

#include <vector>

class CoordTransform {
  typedef double point_t;
  typedef double ref_point_t;
public:
  point_t local_to_global(const ref_point_t& lp, 
      const std::vector<ref_point_t>& lv,
      const std::vector<point_t>& gv) const   {

    point_t gp;
    double lambda[2];
    lambda[0] = (lv[1] - lp)/(lv[1] - lv[0]);
    lambda[1] = (lp - lv[0])/(lv[1] - lv[0]);
    gp = lambda[0]*gv[0] + lambda[1]*gv[1];
    return gp;

  }

  ref_point_t global_to_local(const point_t& gp,
      const std::vector<ref_point_t>& lv,
      const std::vector<point_t>& gv){
    ref_point_t lp;
    double lambda[2];
    lambda[0] = (gv[1] - gp)/(gv[1] - gv[0]);
    lambda[1] = (gp - gv[0])/(gv[1] - gv[0]);
    lp = lambda[0]*lv[0] + lambda[1]*lv[1];
    return lp;
  }
  double local_to_global_jacobian(const ref_point_t& lp,
      const std::vector<ref_point_t>& lv,
      const std::vector<point_t>& gv) const{

    return (gv[1] - gv[0])/(lv[1] - lv[0]);
  }

  double global_to_local_jacobian(const point_t& gp,
      const std::vector<ref_point_t>& lv,
      const std::vector<point_t>& gv) const {

    return (lv[1] - lv[0])/(gv[1] - gv[0]);

  }

};

#endif
