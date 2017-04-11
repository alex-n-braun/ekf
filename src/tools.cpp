#include <iostream>
#include "tools.h"

using namespace std;
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
    VectorXd result(4);
    result<<0,0,0,0;
    {
        auto it_e(estimations.cbegin());
        auto it_g(ground_truth.cbegin());
        for (; it_e!=estimations.cend() & it_g!=ground_truth.cend(); ++it_e, ++it_g) {
            VectorXd diff(*it_e-*it_g);
            diff = diff.array()*diff.array();
            result = result + diff;
        }
    }
    result /= estimations.size();
    result = result.array().sqrt();
    return result;
}

VectorXd Tools::h(const Eigen::VectorXd & x_state) {
    VectorXd result(3);
    float rho(sqrt(x_state[0]*x_state[0]+x_state[1]*x_state[1]));
    float phi;
    if (rho<0.01) {
        rho=0.01;
        phi=0;
    } else
        phi=atan2(x_state[1], x_state[0]);
    float rho_dot((x_state[0]*x_state[2]+x_state[1]*x_state[3])/rho);
    result<<rho, phi, rho_dot;
    return result;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    //check division by zero
    if(fabs(c1) < 0.0001){
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        //return Hj;
        c1 = 0.0001;
    }
    float c2 = sqrt(c1);
    float c3 = (c1*c2);


    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
          -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;
}


