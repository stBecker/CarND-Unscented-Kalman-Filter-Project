#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    cout << "Invalid inputs for RMSE" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for (int i = 0; i < estimations.size(); ++i) {
    // ... your code here
    VectorXd err = estimations[i] - ground_truth[i];
    VectorXd err2 = err.array()*err.array();
    rmse += err2;

  }
  rmse /= estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}