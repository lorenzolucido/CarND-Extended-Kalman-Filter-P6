#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
	P_ = (F_ * P_ * F_.transpose()) + Q_;
}

void KalmanFilter::Update(const VectorXd &z, bool extended) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred;
  VectorXd y;
  MatrixXd H;

  if(extended)
  {
    
    // If both px and py are equal to zero we do not update the Kalman Filter at all.
    // This is a conservative way to avoid division by zero.
    if( x_(0) == 0. and x_(1) == 0.)
      return;

    z_pred = Tools::CartesianToPolar(x_);
    H = Tools::CalculateJacobian(x_); //  Use Jacobian Matrix if non linear
    y = z - z_pred;
    float pi = acos(-1);
    
    if(y[1] < -pi)
      y[1] += 2*pi; 
    if(y[1] > pi)
      y[1] -= 2*pi; 
  }
  else
  {
    H = H_; // Use the standard matrix H_ if linear
    z_pred = H * x_;
    y = z - z_pred;
  }
	MatrixXd Ht = H.transpose();
	MatrixXd S = (H * P_ * Ht) + R_;
	MatrixXd K = P_ * Ht * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
	P_ = (I - K * H) * P_;
}

