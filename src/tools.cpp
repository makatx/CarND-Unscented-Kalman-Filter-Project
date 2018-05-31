#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse = VectorXd(4);
  if (estimations.size() == 0 || ground_truth.size() == 0 || estimations.size()!=ground_truth.size()) {
    cout << "RMSE: Invalid arrays passed" << endl;
    rmse << -1,-1,-1,-1;
    return rmse;
  }
  rmse << 0,0,0,0;
  for(int i=0; i < estimations.size(); ++i){
        VectorXd e = estimations[i];
        VectorXd g = ground_truth[i];
        VectorXd diff = e-g;
        diff = diff.array() * diff.array();
        rmse = diff + rmse;


	}
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();

	return rmse;

}
