#ifndef UKF_H
#define UKF_H

#include <iostream>
#include "Eigen/Dense"
#include "measurement_package.h"

static const float epsilon = 1.0e-05f;
class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // predicted state vector sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_{5};

  // Augmented state dimension
  int n_aug_{7};

  // Sigma point spreading parameter
  double lambda_;

  int64_t previous_timestamp_;  

  public:

  void GenerateAugementedSigmaPoints();
  void SigmaPointPrediction(const double); //pass sigma points through non-linear function
  void PredictMeanAndCovariance();
  bool static floatCompare(float f1, float f2);

  Eigen::MatrixXd Xsig_aug_;
  Eigen::VectorXd x_aug_;
  Eigen::MatrixXd P_aug_;
  Eigen::MatrixXd A_; //for sqrt using L * L.T
  // double px_, py_, v_, yaw_, yawd_, nu_a_, nu_yawdd_;

  
};

#endif  // UKF_H
// Ref. Section
// Formulas used to calculate sigma points and then re-estimating Gausssian from those sigma points can be found under
// http://ais.informatik.uni-freiburg.de/teaching/current-ws/mapping/pdf/slam05-ukf-4.pdf

// Normalized Innovation Squared (NIS) formula used can be found under
// https://web.stanford.edu/group/arl/sites/default/files/public/publications/Robust_TRN_Framework.pdf