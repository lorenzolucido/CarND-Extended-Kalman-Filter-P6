#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  noise_ax = 9.; 
  noise_ay = 9.;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    previous_timestamp_ = measurement_pack.timestamp_;
    VectorXd init_state = VectorXd(4);

    // State covariance matrix
    MatrixXd P_in = MatrixXd(4,4);
    P_in << 1, 0, 0,    0,
			      0, 1, 0,    0,
			      0, 0, 1000, 0,
			      0, 0, 0,    1000;
    
    // State transition matrix (depends on elapsed time - Predict step)
    MatrixXd F_in = MatrixXd(4,4);
    F_in << 1, 0, 1, 0,
			      0, 1, 0, 1,
			      0, 0, 1, 0,
			      0, 0, 0, 1;

    // Process Covariance Matrix (depends on elapsed time  - Predict step)
    MatrixXd Q_in = MatrixXd(4,4);
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      init_state = Tools::PolarToCartesian(measurement_pack.raw_measurements_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      init_state << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      
    }

    // Done initializing, no need to predict or update
    ekf_.Init(init_state, P_in, F_in, H_laser_, R_laser_, Q_in);
    is_initialized_ = true;

    return;
  }

  

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
  float dt2 = pow(dt, 2);
  float dt3 = pow(dt, 3);
  float dt4 = pow(dt, 4);
	ekf_.Q_ << dt4*noise_ax/4,        0,                    dt3*noise_ax/2,       0,
             0,                     dt4*noise_ay/4,       0,                    dt3*noise_ay/2,
	           dt3*noise_ax/2,        0,                    dt2*noise_ax,         0,
	           0,                     dt3*noise_ay/2,       0,                    dt2*noise_ay;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.Update(measurement_pack.raw_measurements_, true);
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_, false);
  }  
}
