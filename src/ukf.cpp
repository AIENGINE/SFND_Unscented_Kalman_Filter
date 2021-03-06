#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial state covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  is_initialized_ = false;
  x_aug_ = VectorXd(n_aug_);
  Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_+1);
  P_aug_ = MatrixXd(n_aug_, n_aug_); //Cholesky (L * L.T) decomposition is performed then sigma points are calculated sigma points (XSig_aug_) using given formulas see ukf.h ref. section
  P_aug_.fill(0.0);

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.0);
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  for (size_t vec_idx = 1; vec_idx < 2*n_aug_+1; ++vec_idx)
  {
    weights_(vec_idx) = 0.5 / (n_aug_ + lambda_);
  }

}

bool UKF::floatCompare(float f1, float f2) {
    if (std::fabs(f1 - f2) <= epsilon)
        return true;
    return std::fabs(f1 - f2) <= epsilon * std::fmax(std::fabs(f1), std::fabs(f2));
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_)
  {

    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0,
            0,
            0;
      P_ << (std_laspx_*std_laspx_), 0, 0, 0, 0,
            0, (std_laspy_*std_laspy_), 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

    }
    else 
    {
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      float px = rho * cos(phi);
      float py = rho * sin(phi); 
      x_ << px,
            py,
            rho_dot,
            phi,
            0;

      P_ << (std_laspx_*std_laspx_), 0, 0, 0, 0,
            0, (std_laspy_*std_laspy_), 0, 0, 0,
            0, 0, (std_radrd_*std_radrd_), 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    }

    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    std::cout<< "..............UKF Initialized.............."<< std::endl;
    std::cout<<"x = " <<std::endl<< x_<<std::endl;
    std::cout<<"P = " <<std::endl<< P_<<std::endl;
  }
  

  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true)
  {
    std::cout<<"Lidar measurement" <<std::endl;
    std::cout<<"x = " <<std::endl<< meas_package.raw_measurements_[0] << " "<< meas_package.raw_measurements_[1] <<std::endl;
    UpdateLidar(meas_package);


  }
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_== true)
  {

    std::cout<<"Radar measurement" <<std::endl;
    UpdateRadar(meas_package);

  }
  else
  {
    std::cout<<"Not a valid measurement";
  }


}

void UKF::Prediction(double delta_t) {

  GenerateAugementedSigmaPoints();
  SigmaPointPrediction(delta_t);
  PredictMeanAndCovariance();

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {

  const int n_z_lidar{2};
  MatrixXd Zsig_lidar = MatrixXd::Zero(n_z_lidar, 2 * n_aug_ + 1);

  VectorXd z_pred_lidar = VectorXd::Zero(n_z_lidar);
  MatrixXd S_lidar = MatrixXd::Zero(n_z_lidar, n_z_lidar);

  MatrixXd R_lidar = MatrixXd(n_z_lidar, n_z_lidar);
  R_lidar << std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;

  VectorXd z_lidar = meas_package.raw_measurements_; //actual measurement from measurement package
  
  MatrixXd T_lidar = MatrixXd::Zero(n_x_, n_z_lidar); // cross correlation matrix

  VectorXd state_diff_vector_lidar = VectorXd::Zero(n_x_); //for calculate cross correlation matrix

  VectorXd measurement_diff_vector_lidar = VectorXd::Zero(n_z_lidar);

  for (int col_idx = 0; col_idx < 2 * n_aug_ + 1; ++col_idx)
  {
      // extract value
    double px = Xsig_pred_(0, col_idx);
    double py = Xsig_pred_(1, col_idx);

    // measurement model
    Zsig_lidar(0, col_idx) = px;   // p_x
    Zsig_lidar(1, col_idx) = py;   // p_y
  }

  // measurement mean state vector
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
      z_pred_lidar += weights_(i) * Zsig_lidar.col(i);
  }

  // System/innovation covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {

    measurement_diff_vector_lidar = Zsig_lidar.col(i) - z_pred_lidar;

    S_lidar += (weights_(i) * measurement_diff_vector_lidar * measurement_diff_vector_lidar.transpose());

    state_diff_vector_lidar = Xsig_pred_.col(i) - x_;

    while (state_diff_vector_lidar(3) > M_PI) state_diff_vector_lidar(3) -= 2.0*M_PI;
    while (state_diff_vector_lidar(3) < -M_PI) state_diff_vector_lidar(3) += 2.0*M_PI;

    T_lidar += (weights_(i) * state_diff_vector_lidar * measurement_diff_vector_lidar.transpose());

  }

  // add measurement noise to system covariance matrix
  S_lidar = S_lidar + R_lidar;


  // for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  // {
  //   // state difference
  //   state_diff_vector_lidar = Xsig_pred_.col(i) - x_;

  //   while (state_diff_vector_lidar(3) > M_PI) state_diff_vector_lidar(3) -= 2.0*M_PI;
  //   while (state_diff_vector_lidar(3) < -M_PI) state_diff_vector_lidar(3) += 2.0*M_PI;

  //   // measurement difference
  //   measurement_diff_vector_lidar = Zsig_lidar.col(i) - z_pred_lidar;

  //   while (measurement_diff_vector_lidar(1) > M_PI) measurement_diff_vector_lidar(1) -= 2.0*M_PI;
  //   while (measurement_diff_vector_lidar(1) < -M_PI) measurement_diff_vector_lidar(1) += 2.0*M_PI;

  //   T_lidar += weights_(i) * state_diff_vector_lidar * measurement_diff_vector_lidar.transpose();
  // }


    // calculate kalman gain K for lidar
    MatrixXd K_lidar = T_lidar * S_lidar.inverse();

    // residuals
    VectorXd z_diff_lidar = z_lidar - z_pred_lidar;

    // angle normalization
    while (z_diff_lidar(1) > M_PI) z_diff_lidar(1) -= 2.0*M_PI;
    while (z_diff_lidar(1) < -M_PI) z_diff_lidar(1) += 2.0*M_PI;

    // update state mean and covariance
    x_ = x_ + K_lidar * z_diff_lidar;
    P_ = P_ - K_lidar * S_lidar * K_lidar.transpose();

  // To test UKF lidar or radar unbaisedness NIS is used here and NIS formula please see Ref. Section ukf.h
  // auto NIS = z_diff_lidar.transpose() * S_lidar.inverse() * z_diff_lidar;
  // std::cout<<"NIS = " <<std::endl<< NIS<<std::endl; 
    

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

  const int n_z_radar{3};
  MatrixXd Zsig_radar = MatrixXd(n_z_radar, 2*n_aug_+1);

  VectorXd z_pred_radar = VectorXd(n_z_radar);
  z_pred_radar.fill(0.0);

  //System/innovation Covariance matrix
  MatrixXd S_radar;
  // Measurement Covariance matrix
  MatrixXd R_radar;

  R_radar = MatrixXd(n_z_radar, n_z_radar);

  R_radar << (std_radr_ * std_radr_), 0, 0,
        0, (std_radphi_ * std_radphi_), 0,
        0, 0, (std_radrd_ * std_radrd_);

  S_radar = MatrixXd::Zero(n_z_radar, n_z_radar);
  // S_radar.fill(0.0);

  // Cross Correlation matrix between predicted sigma points matrix of state and predicted sigma points matrix of measurements
  MatrixXd T_radar = MatrixXd::Zero(n_x_, n_z_radar);
  // T_radar.fill(0.0);
  //actual measurement from measurement package
  VectorXd z_radar = VectorXd::Zero(n_z_radar);
  // z_radar.fill(0.0);

  VectorXd state_diff_vector_radar = VectorXd::Zero(n_x_); //for calculate cross correlation matrix

  VectorXd measurement_diff_vector_radar = VectorXd::Zero(n_z_radar);


  for (int col_idx = 0; col_idx < 2 * n_aug_ + 1; ++col_idx)
  {

    double px = Xsig_pred_(0, col_idx);
    double py = Xsig_pred_(1, col_idx);
    double v  = Xsig_pred_(2, col_idx);
    double yaw = Xsig_pred_(3, col_idx);
    // velocity components for rho_dot_
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig_radar(0, col_idx) = sqrt((px * px) + (py * py));                        // rho
    Zsig_radar(1, col_idx) = atan2(py, px);                                   // phi
    Zsig_radar(2, col_idx) = (px * v1 + py * v2) / sqrt((px * px) + (py * py));   // rho_dot
  }

  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
    z_pred_radar = z_pred_radar + weights_(i) * Zsig_radar.col(i);
  }

  for (int idx = 0; idx < 2 * n_aug_ + 1; ++idx)
  {
    // residual
    measurement_diff_vector_radar = Zsig_radar.col(idx) - z_pred_radar;
    if (measurement_diff_vector_radar(1) > M_PI)
        measurement_diff_vector_radar(1) -= 2 * M_PI;
    if (measurement_diff_vector_radar(1) < -M_PI)
        measurement_diff_vector_radar(1) += 2 * M_PI;
    S_radar = S_radar + (weights_(idx) * measurement_diff_vector_radar * measurement_diff_vector_radar.transpose());

    state_diff_vector_radar = Xsig_pred_.col(idx) - x_;
    if (state_diff_vector_radar(3) > M_PI)
        state_diff_vector_radar(3) -= 2 * M_PI;
    if (state_diff_vector_radar(3) < -M_PI)
        state_diff_vector_radar(3) += 2 * M_PI;

    T_radar += (weights_(idx) * state_diff_vector_radar * measurement_diff_vector_radar.transpose());
  }

  S_radar = S_radar + R_radar;



  z_radar << meas_package.raw_measurements_(0), //rho
        meas_package.raw_measurements_(1),      //phi
        meas_package.raw_measurements_(2);      //rho_dot
  // z_ = meas_package.raw_measurements_;

  // for (size_t col_idx=0; col_idx < 2*n_aug_+1; ++col_idx)
  // {
  //   measurement_diff_vector_ = Zsig_.col(col_idx) - z_pred_;
  //   if (measurement_diff_vector_(1) > M_PI)
  //       measurement_diff_vector_(1) -= 2 * M_PI;
  //   if (measurement_diff_vector_(1) < -M_PI)
  //       measurement_diff_vector_(1) += 2 * M_PI;

  //   state_diff_vector_ = Xsig_pred_.col(col_idx) - x_;
  //   if (state_diff_vector_(3) > M_PI)
  //       state_diff_vector_(3) -= 2 * M_PI;
  //   if (state_diff_vector_(3) < -M_PI)
  //       state_diff_vector_(3) += 2 * M_PI;

  //   T_ += (weights_(col_idx) * state_diff_vector_ * measurement_diff_vector_.transpose());
  // }

  // calculate kalman gain for radar
  MatrixXd K_radar;
  K_radar = T_radar * S_radar.inverse();

    // residual
  VectorXd z_diff_radar = z_radar - z_pred_radar;

  // angle normalization
  while (z_diff_radar(1)> M_PI) z_diff_radar(1)-=2.*M_PI;
  while (z_diff_radar(1)<-M_PI) z_diff_radar(1)+=2.*M_PI;
  x_ = x_ + K_radar * z_diff_radar;
  P_ = P_ - (K_radar * S_radar * K_radar.transpose());

  std::cout<<"x_ = " <<std::endl<< x_<<std::endl;
  std::cout<<"P_ = " <<std::endl<< P_<<std::endl;

  // To test UKF lidar or radar unbaisedness NIS is used here and NIS formula please see Ref. Section ukf.h
  // output values from NIS can then be used to perform consistency check using Chi-squared
  auto NIS = z_diff_radar.transpose() * S_radar.inverse() * z_diff_radar;
  std::cout<<"NIS = " <<std::endl<< NIS<<std::endl; 
}


void UKF::GenerateAugementedSigmaPoints(){

  x_aug_.head(5) = x_;
  x_aug_(5) = 0; 
  x_aug_(6) = 0; 

  P_aug_.topLeftCorner(5, 5) = P_;
  // Augment Q in P_aug_
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6, 6) = std_yawdd_ * std_yawdd_;

  A_ = P_aug_.llt().matrixL();

  // set first column of sigma point matrix
  Xsig_aug_.col(0) = x_aug_;

  // set remaining sigma points
  for (size_t i = 0; i < n_aug_; ++i) {
      Xsig_aug_.col(i+1)     = x_aug_ + sqrt(lambda_ + n_aug_) * A_.col(i); //n_aug_ based idx on the + weighting factor will go from idx 1 to 7
      Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * A_.col(i); //n_aug_ based idx on the - weighting factor will go from idx 8 to 14
  }


}

void UKF::SigmaPointPrediction(const double delta_t){

  // Calculate Predicted state sigma point matrix by passing Xsig_aug_ through non-linear eqs of motion model(CTRV)
  // Xsig_pred represents non-linear state transition from K to K+1

  double dt = delta_t;
  double dt2 = delta_t * delta_t;

  for (size_t col_idx=0; col_idx < 2*n_aug_+1; ++col_idx )
  {

    double px = Xsig_aug_(0, col_idx);
    double py = Xsig_aug_(1, col_idx);
    double v = Xsig_aug_(2, col_idx);
    double yaw = Xsig_aug_(3, col_idx);
    double yawd = Xsig_aug_(4, col_idx);
    double nu_a = Xsig_aug_(5, col_idx);
    double nu_yawdd = Xsig_aug_(6, col_idx);

    if (floatCompare(yawd, 0))
    {
      px = px + v * dt * cos(yaw);
      py = py + v * dt * sin(yaw);
    }
    else
    {
      px = px + (v/yawd) * (sin(yaw + yawd*dt) - sin(yaw));
      py = py + (v/yawd) * (-cos(yaw + yawd*dt) + cos(yaw));
    }
    
    px = px + (0.5 * dt2 * nu_a * cos(yaw));
    py = py + (0.5 * dt2 * nu_a * sin(yaw));

    v = v + (dt * nu_a);
    yaw = yaw + (yawd * dt) + (0.5 * dt2 * nu_yawdd);
    yawd = yawd + (dt * nu_yawdd); 


    Xsig_pred_(0, col_idx) = px;
    Xsig_pred_(1, col_idx) = py;
    Xsig_pred_(2, col_idx) = v;
    Xsig_pred_(3, col_idx) = yaw;
    Xsig_pred_(4, col_idx) = yawd;
    // std::cout<<"Xsig_pred = " <<std::endl<< Xsig_pred_<<std::endl;

  }

}

void UKF::PredictMeanAndCovariance(){

  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);

  VectorXd state_diff_vector = VectorXd::Zero(n_x_);
  
  for (size_t vec_idx = 0; vec_idx < 2*n_aug_+1; ++vec_idx)
  {
    x = x + (weights_(vec_idx) * Xsig_pred_.col(vec_idx));
  }

  for(size_t col_idx=0; col_idx < 2*n_aug_+1; ++col_idx)
  {
    state_diff_vector = Xsig_pred_.col(col_idx) - x;
    if (state_diff_vector(3) > M_PI)
        state_diff_vector(3) -= 2 * M_PI;
    if (state_diff_vector(3) < -M_PI)
        state_diff_vector(3) += 2 * M_PI;
    P = P + (weights_(col_idx) * state_diff_vector * state_diff_vector.transpose());
  }
  x_ = x;
  P_ = P;
  // std::cout<<"x_ = " <<std::endl<< x_<<std::endl;
  // std::cout<<"P_ = " <<std::endl<< P_<<std::endl;

}
