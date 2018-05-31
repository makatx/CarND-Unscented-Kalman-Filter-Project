#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  cout << "In constructor";
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  is_initialized_ = false;
  Xsig_pred_ = MatrixXd(5,2*n_aug_+1);
  weights_ = VectorXd(2*n_aug_+1);

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  cout << "In ProcessMeasurment..."<<endl;
  if(!is_initialized_) {
    cout << "Initializing State..."<<endl;
    time_us_ = meas_package.timestamp_;
    float px, py, v, phi, phi_dot;
    P_.fill(0.0);
    if(meas_package.sensor_type_==MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);
      px = rho * cos(theta);
      py = rho * sin(theta);
      /*although RADAR does not measure the object's actual velocity,
      it does however measure the radial component the same which is
      always a fraction (cos(theta-phi)) of the actual velocity.
      since we dont have other info at initialization, we'll set the velocity of state vector
      to rho_dot */
      v = rho_dot;
      phi = 0;
      phi_dot=0;
      P_(0,0) = 0.027;
      P_(1,1) = 0.027;
      P_(2,2) = 10;
      P_(3,3) = 10;
      P_(4,4) = 10;

    }
    else {
      px = meas_package.raw_measurements_(0);
      py = meas_package.raw_measurements_(1);
      v = 0;
      phi = 0;
      phi_dot=0;
      P_(0,0) = std_laspx_*std_laspx_;
      P_(1,1) = std_laspy_*std_laspy_;
      P_(2,2) = 10;
      P_(3,3) = 10;
      P_(4,4) = 10;
    }
    x_ << px, py, v, phi, phi_dot;
    is_initialized_ = true;
    cout << "State initialized"<<endl;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  cout << "About to predict..."<<endl;
  Prediction(delta_t);

  cout << "About to update..."<<endl;
  if(meas_package.sensor_type_==MeasurementPackage::RADAR && use_radar_)
    UpdateRadar(meas_package);
  else if(meas_package.sensor_type_==MeasurementPackage::LASER && use_laser_)
    UpdateLidar(meas_package);

  return;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  AugmentedSigmaPoints(&Xsig_aug);
  SigmaPointPrediction(&Xsig_pred_, Xsig_aug, delta_t);
  PredictMeanAndCovariance(&x_, &P_);

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 3;
  MatrixXd H_ = MatrixXd(2,5);
  H_ << 1,0,0,0,0,
        0,1,0,0,0;
  MatrixXd R_ = MatrixXd(2,2);
  R_ << std_laspx_*std_laspx_,0,
        0, std_laspy_*std_laspy_;

  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  VectorXd y = z - H_*x_;
  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ = x_ + K*y;
  MatrixXd I = MatrixXd::Identity(5,5);
  P_ = (I - K*H_)*P_;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //create matrix for sigma points in measurement space
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //Update RADAR
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();


}

/**
 * AugmentedSigmaPoints
 * @param Xsig_out The variable to place generated augumented sigma points in
 * @param delta_t seconds past since last call
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  *Xsig_out = Xsig_aug;

}

/**
 * SigmaPointPrediction
 * @param Xsig_out The variable to place prediction for sigma points in
 * @param Xsig_aug Augumented sigma point matrix
 * @param delta_t seconds past since last call
 */
void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, MatrixXd Xsig_aug, double delta_t) {

  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  //predict sigma points
  for(int i = 0; i < 2*n_aug_+1; i++) {
      double px = Xsig_aug(0,i);
      double py = Xsig_aug(1, i);
      double v = Xsig_aug(2, i);
      double phi = Xsig_aug(3, i);
      double phi_dot = Xsig_aug(4, i);
      double va = Xsig_aug(5, i);
      double vp = Xsig_aug(6, i);

      if(fabs(phi_dot)>0.001) {
          Xsig_pred(0,i) = px + (v/phi_dot) * (sin(phi+phi_dot*delta_t)-sin(phi)) + 0.5*delta_t*delta_t*cos(phi)*va;
          Xsig_pred(1,i) = py + (v/phi_dot) * (-cos(phi+phi_dot*delta_t)+cos(phi)) + 0.5*delta_t*delta_t*sin(phi)*va;
          Xsig_pred(2,i) = v + delta_t*va;
          Xsig_pred(3,i) = phi + phi_dot*delta_t + 0.5*delta_t*delta_t*vp;
          Xsig_pred(4,i) = phi_dot + delta_t*vp;
      }
      else {
          Xsig_pred(0,i) = px + v*cos(phi)*delta_t + 0.5*delta_t*delta_t*cos(phi)*va;
          Xsig_pred(1,i) = py + v*sin(phi)*delta_t + 0.5*delta_t*delta_t*sin(phi)*va;
          Xsig_pred(2,i) = v + delta_t*va;
          Xsig_pred(3,i) = phi + phi_dot*delta_t + 0.5*delta_t*delta_t*vp;
          Xsig_pred(4,i) = phi_dot + delta_t*vp;
      }
    }
    *Xsig_out = Xsig_pred;

}

/**
 * PredictMeanAndCovariance
 * @param x_pred Mean of predicted points
 * @param P_pred covarianceof predicted points
 */
void UKF::PredictMeanAndCovariance(VectorXd* x_pred, MatrixXd* P_pred) {

  //create vector for predicted state
 VectorXd x = VectorXd(n_x_);

 //create covariance matrix for prediction
 MatrixXd P = MatrixXd(n_x_, n_x_);

 //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  *x_pred = x;
  *P_pred = P;

}
