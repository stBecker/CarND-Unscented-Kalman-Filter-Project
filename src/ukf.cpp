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
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.5;
  
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
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  P_ = MatrixXd::Identity(n_x_, n_x_);

  //calculate innovation covariance matrix S
  R_radar = MatrixXd::Zero(3, 3);
  R_radar(0, 0) = std_radr_*std_radr_;
  R_radar(1, 1) = std_radphi_*std_radphi_;
  R_radar(2, 2) = std_radrd_*std_radrd_;

  R_laser = MatrixXd::Zero(2, 2);
  R_laser(0, 0) = std_laspx_*std_laspx_;
  R_laser(1, 1) = std_laspy_*std_laspy_;

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
  if (!is_initialized_) {
    cout << "Initializing" << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float ro_dot;
      float x = ro * cosf(theta);
      float y = ro * sinf(theta);
      x_ << x, y, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (!use_radar_) {
      return;
    }
  }
  else {
    if (!use_laser_) {
      return;
    }
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else {
    UpdateLidar(meas_package);
  }

  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
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
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd xk = Xsig_aug.col(i).head(n_x_);
    float v = xk(2);
    float psi = xk(3);
    float psi_dot = xk(4);
    float nu_a = Xsig_aug(n_x_, i);
    float nu_psi_dotdot = Xsig_aug(n_x_ + 1, i);

    VectorXd xk_dot = VectorXd(n_x_);
    //avoid division by zero
    if (fabs(psi_dot) < 0.0001) {
      xk_dot <<
        v*cos(psi)*delta_t,
        v*sin(psi)*delta_t,
        0,
        0,
        0;
    }
    else {
      xk_dot <<
        v / psi_dot*(sin(psi + psi_dot*delta_t) - sin(psi)),
        v / psi_dot*(-cos(psi + psi_dot*delta_t) + cos(psi)),
        0,
        psi_dot*delta_t,
        0;
    }

    VectorXd xk_nu = VectorXd(n_x_);
    xk_nu <<
      0.5*delta_t*delta_t*cos(psi)*  nu_a,
      0.5*delta_t*delta_t*sin(psi)*  nu_a,
      delta_t*nu_a,
      0.5*delta_t*delta_t* nu_psi_dotdot,
      delta_t* nu_psi_dotdot;

    VectorXd xk1 = xk + xk_dot + xk_nu;

    //write predicted sigma points into right column
    Xsig_pred_.col(i) = xk1;
  }

  //predict state mean
  x_ = Xsig_pred_ * weights_;

  //predict state covariance matrix
  P_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    MatrixXd r = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (r(3)> M_PI) r(3) -= 2.*M_PI;
    while (r(3)<-M_PI) r(3) += 2.*M_PI;

    MatrixXd rt = r.transpose();
    MatrixXd p = weights_(i)*r*rt;
    P_ += p;
  }
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
  int n_z = 2;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  ////create matrix for sigma points in measurement space
  //MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  ////transform sigma points into measurement space
  //for (int i = 0; i < 2 * n_aug_ + 1; i++) {
  //  float px = Xsig_pred_(0, i);
  //  float py = Xsig_pred_(1, i);

  //  Zsig.col(i) << px, py;
  //}

  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, 2 * n_aug_ + 1);

  //calculate mean predicted measurement
  z_pred = Zsig*weights_;

  S.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    MatrixXd r = Zsig.col(i) - z_pred;
    S += weights_(i)*r*r.transpose();
  }
  S += R_laser;

  VectorXd z = meas_package.raw_measurements_;
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc * Si;

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();

  double nis = z_diff.transpose()*Si*z_diff;
  cout << "NIS: " << nis << endl;

  if (nis >= confidence_2df) {
    count_laser_outlier++;
  }
  count_laser++;

  cout << "Laser outliers: " << ((double)count_laser_outlier) / count_laser << endl;
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
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    float px = Xsig_pred_(0, i);
    float py = Xsig_pred_(1, i);
    float v = Xsig_pred_(2, i);
    float sci = Xsig_pred_(3, i);
    float scid = Xsig_pred_(4, i);

    if (px == 0 && py == 0) {
      px = 0.0001;
    }
    float rho = sqrt(px*px + py*py);
    float phi = atan2(py, px);
    float rhod = (px*cos(sci)*v + py*sin(sci)*v) / rho;
    Zsig.col(i) << rho, phi, rhod;
  }

  //calculate mean predicted measurement
  z_pred = Zsig*weights_;

  S.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    MatrixXd r = Zsig.col(i) - z_pred;
    //angle normalization
    while (r(1)> M_PI) r(1) -= 2.*M_PI;
    while (r(1)<-M_PI) r(1) += 2.*M_PI;
    S += weights_(i)*r*r.transpose();
  }
  S += R_radar;

  VectorXd z = meas_package.raw_measurements_;
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc * Si;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();

  double nis = z_diff.transpose()*Si*z_diff;
  cout << "NIS: " << nis << endl;

  if (nis >= confidence_3df) {
    count_radar_outlier++;
  }
  count_radar++;

  cout << "Radar outliers: " << ((double)count_radar_outlier)/count_radar << endl;

}
