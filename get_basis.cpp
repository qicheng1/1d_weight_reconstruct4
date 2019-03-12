#include <vector>
#include "Eigen/Dense"
std::vector<std::vector<double>> get_basis(Eigen::MatrixXd& T, Eigen::MatrixXd& Tt)
{
    int i,j,k;
    std::vector<std::vector<double>> mtx(Tt.cols(),std::vector<double>(Tt.rows()));
    for(i = 0; i<T.rows(); i++)
	for(j = 0; j<Tt.cols(); j++)
	    for(k = 0; k<T.cols(); k++)
		mtx[j][i] = mtx[j][i]+T(i,k)*Tt(k,j);
    return mtx;
}

std::vector<std::vector<double>> get_basis(Eigen::MatrixXd& T,
		Eigen::MatrixXd& Tt,
		Eigen::MatrixXd& W)
{
    int i,j,k;
    std::vector<std::vector<double>> mtx(Tt.cols(),std::vector<double>(Tt.rows()));
    Eigen::MatrixXd temp = W.transpose()*Tt.transpose();
    for(i = 0; i<temp.rows(); i++)
        for(j = 0; j<T.cols(); j++)
            for(k = 0; k<temp.cols(); k++)
                mtx[i][j] = mtx[i][j]+temp(i,k)*T(j,k);
    return mtx;

}
