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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  is_initialized_ = false;
  x_aug_ = VectorXd(n_aug_);
  Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_+1);
  P_aug_ = MatrixXd(n_aug_, n_aug_); //Cholesky (L * L.T) decomposition is performed then sigma points are calculated sigma points (XSig_aug_) using given formulas see ukf.h ref. section
  // A_ = MatrixXd(n_aug_, n_aug_);
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.0);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  float weights_value = 0.5 / (n_aug_ + lambda_);
  for (size_t vec_idx = 1; vec_idx < 2*n_aug_+1; ++vec_idx)
  {
    weights_(vec_idx) = 0.5 / (n_aug_ + lambda_);
  }

  px_ = py_ = v_ = yaw_ = yawd_= nu_a_ = nu_yawdd_ = 0.0;

}

bool UKF::floatCompare(float f1, float f2) {
    if (std::fabs(f1 - f2) <= epsilon)
        return true;
    return std::fabs(f1 - f2) <= epsilon * std::fmax(std::fabs(f1), std::fabs(f2));
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: based on incoming measurement from meas_package which has lidar and radar initizalized the x_ and P_
   * Calcualte the dt from the last and the current timesstamp -> conver the timestamp into seconds
   * 
   */
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
    else // toggle this part in stepHighway function
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
            0, 0, 0, (std_radphi_*std_radphi_), 0,
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

  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    std::cout<<"Lidar measurement";
    std::cout<<"x = " <<std::endl<< meas_package.raw_measurements_[0] << " "<< meas_package.raw_measurements_[1] <<std::endl;
    UpdateLidar(meas_package);


  }
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {

    UpdateRadar(meas_package);

  }
  else
  {
    std::cout<<"Not a valid measurement";
  }


}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Augment and Predict sigma points x_aug, P_aug, Xsig_aug, then Predict Sigma points Xsig_pred, 
   * and then predict state covariance matrix and state vector x_ and P_. This step remain same for lidar and radar as 
   * mean state vector is the same for both of the sensors. 
   */
  GenerateAugementedSigmaPoints();
  SigmaPointPrediction(delta_t);

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS(Normalized Innovation Squared), if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}



void UKF::GenerateAugementedSigmaPoints(){

  x_aug_.head(5) = x_;
  x_aug_(5) = 0; 
  x_aug_(6) = 0; 

  P_aug_.topLeftCorner(5, 5) = P_;
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  A_ = P_aug_.llt().matrixL();

  // set first column of sigma point matrix
  Xsig_aug_.col(0) = x_aug_;

  // set remaining sigma points
  for (int i = 0; i < n_aug_; ++i) {
      Xsig_aug_.col(i+1)     = x_aug_ + sqrt(lambda_ + n_aug_) * A_.col(i); //n_aug_ based idx on the + weighting factor will go from idx 1 to 7
      Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * A_.col(i); //n_aug_ based idx on the - weighting factor will go from idx 8 to 14
  }


}

void UKF::SigmaPointPrediction(const double delta_t){

  // Predict state sigma point matrix by passing Xsig_aug_ throught non-linear eqs of motion model(CTRV)

  double dt = delta_t;
  double dt2 = delta_t * delta_t;

  for (size_t col_idx=0; col_idx < 2*n_aug_+1; ++col_idx )
  {

    px_ = Xsig_aug_(0, col_idx);
    py_ = Xsig_aug_(1, col_idx);
    v_ = Xsig_aug_(2, col_idx);
    yaw_ = Xsig_aug_(3, col_idx);
    yawd_ = Xsig_aug_(4, col_idx);
    nu_a_ = Xsig_aug_(5, col_idx);
    nu_yawdd_ = Xsig_aug_(6, col_idx);

    if (floatCompare(yawd_, 0))
    {
      px_ = px_ + v_ * dt * cos(yaw_);
      py_ = py_ + v_ * dt * sin(yaw_);
    }
    else
    {
      px_ = px_ + (v_/yawd_) * (sin(yaw_ + yawd_*dt) - sin(yaw_));
      py_ = py_ + (v_/yawd_) * (-cos(yaw_ + yawd_*dt) + cos(yaw_));
    }
    
    px_ = px_ + (0.5 * dt2 * nu_a_ * cos(yaw_));
    py_ = py_ + (0.5 * dt2 * nu_a_ * sin(yaw_));

    v_ = v_ + (dt * nu_a_);
    yaw_ = yaw_ + (yawd_ * dt) + (0.5 * dt2 * nu_yawdd_);
    yawd_ = yawd_ + (dt * nu_yawdd_); 


    Xsig_pred_(0,col_idx) = px_;
    Xsig_pred_(1,col_idx) = py_;
    Xsig_pred_(2,col_idx) = v_;
    Xsig_pred_(3,col_idx) = yaw_;
    Xsig_pred_(4,col_idx) = yawd_;
    std::cout<<"Xsig_pred = " <<std::endl<< Xsig_pred_<<std::endl;

  }

}

