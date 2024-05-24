#ifndef HYDRO_SOLVER_METHODS_LAGRANGIAN_FEM_METHOD_HPP_
#define HYDRO_SOLVER_METHODS_LAGRANGIAN_FEM_METHOD_HPP_

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../method.hpp"

class LagrangianFemMethod : public Method {
 public:
  LagrangianFemMethod(const Problem &problem_, const std::size_t xCells_,
                      const std::size_t yCells_, const std::size_t order_);

  LagrangianFemMethod(const LagrangianFemMethod &) = delete;
  LagrangianFemMethod(LagrangianFemMethod &&) = delete;

  LagrangianFemMethod &operator=(const LagrangianFemMethod &) = delete;
  LagrangianFemMethod &operator=(LagrangianFemMethod &&) = delete;

  virtual ~LagrangianFemMethod() = default;

 public:
  virtual void dumpSolverInfo() const override;
  virtual void dumpData() const override;

  virtual void calcStep() override;

 private:
  // STATIC INITIALIZATION FUNCTIONS
  static double initHminCoeff(
      const std::size_t numberOfKinematicPointsPerCellPerDimention);
  static std::vector<double> initQuadratureWeights(
      const std::size_t numberOfQuadPointsPerDimention,
      const std::size_t solverDimention);
  static Eigen::MatrixXd initKinematicBasisQuadValues(
      const std::size_t numberOfKinematicPointsPerCellPerDimention,
      const std::size_t numberOfQuadraturePointsPerCellPerDimention,
      const std::size_t solverDimention);
  static Eigen::MatrixXd initKinematicBasisdxQuadValues(
      const std::size_t numberOfKinematicPointsPerCellPerDimention,
      const std::size_t numberOfQuadraturePointsPerCellPerDimention,
      const std::size_t solverDimention);
  static Eigen::MatrixXd initKinematicBasisdyQuadValues(
      const std::size_t numberOfKinematicPointsPerCellPerDimention,
      const std::size_t numberOfQuadraturePointsPerCellPerDimention,
      const std::size_t solverDimention);
  static Eigen::MatrixXd initThermodynamicBasisQuadValues(
      const std::size_t numberOfThermodynamicPointsPerCellPerDimention,
      const std::size_t numberOfQuadraturePointsPerCellPerDimention,
      const std::size_t solverDimention);
  // INITIALIZATION FUNCTIONS
  void initKinematicVectors();
  void initQuadVectors();
  void initVolFracVector();
  void initRhoVector();
  void initThermodynamicVector();
  void initKinematicMassMatrix();
  void initThermoMassMatrixInv();
  void initForceMatrix();
  void initVolFracRateVector();
  void initKinematicSolver();

  // CALC METHODS
  void calcKinematicMassMatrix();
  void calcThermoMassMatrixInv();
  void calcForceMatrix(const Eigen::VectorXd &xCalc,
                       const Eigen::VectorXd &uCalc,
                       const Eigen::VectorXd &volFracCalc,
                       const Eigen::VectorXd &eCalc);

  // QUAD METHODS
  Eigen::MatrixXd quadKinematicCellMass(const std::size_t cell);
  Eigen::MatrixXd quadThermoCellMass(const std::size_t cell,
                                     const std::size_t material);
  Eigen::MatrixXd quadForceMatrix(const std::size_t cell,
                                  const Eigen::VectorXd &xCalc,
                                  const Eigen::VectorXd &uCalc,
                                  const Eigen::VectorXd &volFracCalc,
                                  const Eigen::VectorXd &eCalc);

  // STRESS CALC METHODS
  Eigen::Matrix2d calcStressTensor(
      const std::size_t cell, const std::size_t quad,
      const std::size_t material, double &soundSpeed, double &rhoLocal,
      double &pLocal, double &maxViscosityCoeff,
      const Eigen::Matrix2d &jacobian, const double jacobianDet,
      const Eigen::Matrix2d &jacobianInv, const Eigen::VectorXd &uCalc,
      const Eigen::VectorXd &volFracCalc, const Eigen::VectorXd &eCalc);
  Eigen::Matrix2d calcArtificialViscosity(
      const std::size_t cell, const std::size_t quad, const double soundSpeed,
      const double rhoLocal, double &maxViscosityCoeff,
      const Eigen::Matrix2d &jacobian, const Eigen::Matrix2d &jacobianInv,
      const Eigen::Matrix2d &jacobianInitial, const Eigen::VectorXd &uCalc);
  double calcViscosityCoeff(const double soundSpeed, const double rhoLocal,
                            const double eigenvalue,
                            const Eigen::Vector2d &eigenvector,
                            const Eigen::Matrix2d &velocityGrad,
                            const Eigen::Matrix2d &jacobian,
                            const Eigen::Matrix2d &jacobianInitial);
  void calcTau(double hmin, double soundSpeed, double rhoLocal,
               double maxViscosityCoeff);

  // BOUNDARY RESOLVE FUNCTIONS
  void resolveBoundaryKinematicMassMatrix();
  void resolveLeftBoundaryKinematicMassMatrix();
  void resolveTopBoundaryKinematicMassMatrix();
  void resolveRightBoundaryKinematicMassMatrix();
  void resolveBottomBoundaryKinematicMassMatrix();

  void resolveBoundaryForceVector(Eigen::VectorXd &FuCalc);
  void resolveLeftBoundaryForceVector(Eigen::VectorXd &FuCalc);
  void resolveTopBoundaryForceVector(Eigen::VectorXd &FuCalc);
  void resolveRightBoundaryForceVector(Eigen::VectorXd &FuCalc);
  void resolveBottomBoundaryForceVector(Eigen::VectorXd &FuCalc);

  // STEP FUNCTIONS
  void RK2Step();

  // UTILITY FUNCTIONS
  std::size_t getKinematicIndexFromCell(const std::size_t cell,
                                        const std::size_t node,
                                        const std::size_t direction) const;
  std::size_t getThermodynamicIndexFromCell(const std::size_t cell,
                                            const std::size_t node,
                                            const std::size_t material) const;
  std::size_t getQuadIndexFromCell(const std::size_t cell,
                                   const std::size_t quad,
                                   const std::size_t material) const;
  std::size_t getKinematicIndexLocal(const std::size_t node,
                                     const std::size_t direction) const;
  std::size_t getThermoIndexLocal(const std::size_t node,
                                  const std::size_t material) const;
  Eigen::Matrix2d calcJacobian(const std::size_t cell, const std::size_t quad,
                               const Eigen::VectorXd &xCalc) const;

 private:
  const std::size_t xCells;
  const std::size_t yCells;
  const std::size_t kNumberOfCells;
  const std::size_t kOrder;
  const std::size_t kSolverDimention = 2;
  const std::size_t kNumberOfMaterials;

  const std::size_t kNumberOfKinematicPointsPerCellPerDimention;
  const std::size_t kNumberOfThermodynamicPointsPerCellPerDimention;
  const std::size_t kNumberOfQuadraturePointsPerCellPerDimention;
  const std::size_t kNumberOfKinematicPointsPerCell;
  const std::size_t kNumberOfThermodynamicPointsPerCell;
  const std::size_t kNumberOfQuadraturePointsPerCell;
  const std::size_t kNumberOfKinematicPointsTotal;
  const std::size_t kNumberOfThermodynamicPointsTotal;
  const std::size_t kNumberOfQuadraturePointsTotal;

  const double q1 = 0.5;
  const double q2 = 2.0;
  const double alpha = 0.5;
  const double alphamu = 2.5;
  const double beta1 = 0.85;
  const double beta2 = 1.02;
  const double gamma = 0.8;
  const double hminCoeff;

  double l0;

  double tau;

  const std::vector<double> quadWeights;
  const Eigen::MatrixXd kinematicBasisQuadValues;
  const Eigen::MatrixXd kinematicBasisdxQuadValues;
  const Eigen::MatrixXd kinematicBasisdyQuadValues;
  const Eigen::MatrixXd thermodynamicBasisQuadValues;

  Eigen::VectorXd x;           // Kinematic
  Eigen::VectorXd u;           // Kinematic
  Eigen::VectorXd volFrac;     // Quad
  Eigen::VectorXd rhoInitial;  // Quad
  Eigen::VectorXd e;           // Thermo

  Eigen::VectorXd xInitial;        // Kinematic
  Eigen::VectorXd volFracInitial;  // Quad

  Eigen::SparseMatrix<double> kinematicMassMatrix;
  Eigen::SparseMatrix<double> thermoMassMatrixInv;

  Eigen::SparseMatrix<double> forceMatrix;
  Eigen::VectorXd volFracRate;

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      kinematicMassMatrixSolver;

  // RK2 VARIABLES
  Eigen::VectorXd x05;
  Eigen::VectorXd u05;
  Eigen::VectorXd volFrac05;
  Eigen::VectorXd e05;

  Eigen::VectorXd Fu;
  Eigen::VectorXd Fe;
};

inline void LagrangianFemMethod::calcTau(double hmin, double soundSpeed,
                                         double rhoLocal,
                                         double maxViscosityCoeff) {
  const double denominator = soundSpeed / hmin + alphamu * maxViscosityCoeff /
                                                     (rhoLocal * hmin * hmin);
  const double tauLocal = alpha / denominator;

  if (tauLocal < tau) {
    tau = tauLocal;
  }
}

inline std::size_t LagrangianFemMethod::getKinematicIndexFromCell(
    const std::size_t cell, const std::size_t node,
    const std::size_t direction) const {
  assert(cell < kNumberOfCells);
  assert(node < kNumberOfKinematicPointsPerCell);
  assert(direction < kSolverDimention);

  const std::size_t celli = cell / yCells;
  const std::size_t cellj = cell % yCells;

  const std::size_t nodei = node / kNumberOfKinematicPointsPerCellPerDimention;
  const std::size_t nodej = node % kNumberOfKinematicPointsPerCellPerDimention;

  return kSolverDimention * ((celli * kOrder + nodei) * (yCells * kOrder + 1) +
                             cellj * kOrder + nodej) +
         direction;
}

inline std::size_t LagrangianFemMethod::getThermodynamicIndexFromCell(
    const std::size_t cell, const std::size_t node,
    const std::size_t material) const {
  assert(cell < kNumberOfCells);
  assert(node < kNumberOfThermodynamicPointsPerCell);
  assert(material < kNumberOfMaterials);

  return kNumberOfMaterials *
             (cell * kNumberOfThermodynamicPointsPerCell + node) +
         material;
}

inline std::size_t LagrangianFemMethod::getQuadIndexFromCell(
    const std::size_t cell, const std::size_t quad,
    const std::size_t material) const {
  assert(cell < kNumberOfCells);
  assert(quad < kNumberOfQuadraturePointsPerCell);
  assert(material < kNumberOfMaterials);

  return kNumberOfMaterials * (cell * kNumberOfQuadraturePointsPerCell + quad) +
         material;
}

inline std::size_t LagrangianFemMethod::getKinematicIndexLocal(
    const std::size_t node, const std::size_t direction) const {
  assert(node < kNumberOfKinematicPointsPerCell);
  assert(direction < kSolverDimention);
  return kSolverDimention * node + direction;
}

inline std::size_t LagrangianFemMethod::getThermoIndexLocal(
    const std::size_t node, const std::size_t material) const {
  assert(node < kNumberOfThermodynamicPointsPerCell);
  assert(material < kNumberOfMaterials);
  return kNumberOfMaterials * node + material;
}

#endif  // HYDRO_SOLVER_METHODS_LAGRANGIAN_FEM_METHOD_HPP_
