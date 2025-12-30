#include "plastic_elastic_fem_method.hpp"

#include <Eigen/src/Core/Matrix.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "../utils/basis.hpp"
#include "../utils/legendre_nodes.hpp"
#include "../utils/legendre_weights.hpp"
#include "../utils/lobatto_nodes.hpp"

PlasticElasticMethod::PlasticElasticMethod(const Problem &problem_,
                                           const std::size_t xCells_,
                                           const std::size_t yCells_,
                                           const std::size_t order_)
    : Method(std::to_string(order_) + "-" + std::to_string(xCells_) + "x" +
                 std::to_string(yCells_),
             problem_),
      xCells(xCells_),
      yCells(yCells_),
      kNumberOfCells(xCells * yCells),
      kOrder(order_),
      kNumberOfMaterials(problem.numberOfMaterials),
      kNumberOfKinematicPointsPerCellPerDimention(kOrder + 1),
      kNumberOfThermodynamicPointsPerCellPerDimention(kOrder),
      kNumberOfQuadraturePointsPerCellPerDimention(2 * kOrder),

      // kNumberOfOutputPointsPerCellPerDimention(kOrder == 1 ? 1
      //                                                      : 5 * (kOrder /
      //                                                      2)),

      kNumberOfOutputPointsPerCellPerDimention(
          kNumberOfQuadraturePointsPerCellPerDimention),

      kNumberOfKinematicPointsPerCell(std::pow(
          kNumberOfKinematicPointsPerCellPerDimention, kSolverDimention)),
      kNumberOfThermodynamicPointsPerCell(std::pow(
          kNumberOfThermodynamicPointsPerCellPerDimention, kSolverDimention)),
      kNumberOfQuadraturePointsPerCell(std::pow(
          kNumberOfQuadraturePointsPerCellPerDimention, kSolverDimention)),
      kNumberOfOutputPointsPerCell(
          std::pow(kNumberOfOutputPointsPerCellPerDimention, kSolverDimention)),
      kNumberOfKinematicPointsTotal((kOrder * xCells + 1) *
                                    (kOrder * yCells + 1)),
      kNumberOfThermodynamicPointsTotal(kNumberOfThermodynamicPointsPerCell *
                                        kNumberOfCells),
      kNumberOfQuadraturePointsTotal(kNumberOfQuadraturePointsPerCell *
                                     kNumberOfCells),
      kNumberOfOutputPointsTotal(kNumberOfOutputPointsPerCell * kNumberOfCells),
      hminCoeff(initHminCoeff(kNumberOfKinematicPointsPerCellPerDimention)),
      quadWeights(initQuadratureWeights(
          kNumberOfQuadraturePointsPerCellPerDimention, kSolverDimention)),
      kinematicBasisQuadValues(initKinematicBasisQuadValues(
          kNumberOfKinematicPointsPerCellPerDimention,
          kNumberOfQuadraturePointsPerCellPerDimention, kSolverDimention)),
      kinematicBasisdxQuadValues(initKinematicBasisdxQuadValues(
          kNumberOfKinematicPointsPerCellPerDimention,
          kNumberOfQuadraturePointsPerCellPerDimention, kSolverDimention)),
      kinematicBasisdyQuadValues(initKinematicBasisdyQuadValues(
          kNumberOfKinematicPointsPerCellPerDimention,
          kNumberOfQuadraturePointsPerCellPerDimention, kSolverDimention)),
      thermodynamicBasisQuadValues(initThermodynamicBasisQuadValues(
          kNumberOfThermodynamicPointsPerCellPerDimention,
          kNumberOfQuadraturePointsPerCellPerDimention, kSolverDimention)),

      // kinematicBasisOutputValues(initKinematicBasisOutputValues(
      //     kNumberOfKinematicPointsPerCellPerDimention,
      //     kNumberOfOutputPointsPerCellPerDimention, kSolverDimention)),
      // kinematicBasisdxOutputValues(initKinematicBasisdxOutputValues(
      //     kNumberOfKinematicPointsPerCellPerDimention,
      //     kNumberOfOutputPointsPerCellPerDimention, kSolverDimention)),
      // kinematicBasisdyOutputValues(initKinematicBasisdyOutputValues(
      //     kNumberOfKinematicPointsPerCellPerDimention,
      //     kNumberOfOutputPointsPerCellPerDimention, kSolverDimention)),
      // thermodynamicBasisOutputValues(initThermodynamicBasisOutputValues(
      //     kNumberOfThermodynamicPointsPerCellPerDimention,
      //     kNumberOfOutputPointsPerCellPerDimention, kSolverDimention)),
      // quadBasisOutputValues(initQuadBasisOutputValues(
      //     kNumberOfQuadraturePointsPerCellPerDimention,
      //     kNumberOfOutputPointsPerCellPerDimention, kSolverDimention))

      kinematicBasisOutputValues(initKinematicBasisQuadValues(
          kNumberOfKinematicPointsPerCellPerDimention,
          kNumberOfOutputPointsPerCellPerDimention, kSolverDimention)),
      kinematicBasisdxOutputValues(initKinematicBasisdxQuadValues(
          kNumberOfKinematicPointsPerCellPerDimention,
          kNumberOfOutputPointsPerCellPerDimention, kSolverDimention)),
      kinematicBasisdyOutputValues(initKinematicBasisdyQuadValues(
          kNumberOfKinematicPointsPerCellPerDimention,
          kNumberOfOutputPointsPerCellPerDimention, kSolverDimention)),
      thermodynamicBasisOutputValues(initThermodynamicBasisQuadValues(
          kNumberOfThermodynamicPointsPerCellPerDimention,
          kNumberOfOutputPointsPerCellPerDimention, kSolverDimention)),
      quadBasisOutputValues(Eigen::MatrixXd::Identity(
          kNumberOfOutputPointsPerCell, kNumberOfOutputPointsPerCell))

{
  initKinematicVectors();
  initQuadVectors();
  initThermodynamicVector();
  initKinematicMassMatrix();
  initThermoMassMatrixInv();
  initForceMatrix();
  initVolFracRateVector();
  initEnergyExchangeRateVector();
  initDeviatorStressTensorVectors();
  initDeviatorStressRateVectors();
  initEquivalentPlasticStrainVector();
  initKinematicSolver();
}

void PlasticElasticMethod::dumpSolverInfo() const {
  {
    const std::string solverInfoFilename = "solver_info.txt";
    std::ofstream ofs(outputDirPath / solverInfoFilename);

    if (!ofs) {
      const std::string error_message =
          "FAILED TO OPEN FILE: " +
          (outputDirPath / solverInfoFilename).string();

      std::cout << error_message << std::endl;
      return;
    }

    ofs << "XCELLS " << std::to_string(xCells) << "\n";
    ofs << "YCELLS " << std::to_string(yCells) << "\n";
    ofs << "ORDER " << std::to_string(kOrder) << "\n";
    ofs << "DIMENSIONS " << std::to_string(kSolverDimention) << "\n";
    ofs << "MATERIALS " << std::to_string(kNumberOfMaterials) << "\n";
    ofs << "IMAX " << kNumberOfOutputPointsPerCellPerDimention * xCells << "\n";
    ofs << "JMAX " << kNumberOfOutputPointsPerCellPerDimention * yCells << "\n";
    ofs << "OUT_PER_CELL_PER_DIM " << kNumberOfOutputPointsPerCellPerDimention
        << "\n";
    ofs << "OUT_PER_CELL " << kNumberOfOutputPointsPerCell << "\n";
  }
}

void PlasticElasticMethod::dumpData() const {
  const std::size_t tmp = std::round(1000.0 * t * problem.tMul);
  const std::size_t integralPart = tmp / 1000;
  const std::size_t decimalPart = tmp % 1000;
  const std::filesystem::path dataOutputDir =
      outputDirPath /
      (std::to_string(integralPart) + "_" + std::to_string(decimalPart));

  std::filesystem::create_directory(dataOutputDir);

  // DUMP X
  {
    std::ofstream ofs(dataOutputDir / "x.fem", std::ios::binary);
    constexpr std::size_t direction = 0;
    for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
      for (std::size_t out = 0; out < kNumberOfOutputPointsPerCell; out++) {
        double xOut = 0.0;
        for (std::size_t node = 0; node < kNumberOfKinematicPointsPerCell;
             node++) {
          const std::size_t nodeIndex =
              getKinematicIndexFromCell(cell, node, direction);
          const double xLocal = x(nodeIndex);
          const double nodeBasis = kinematicBasisOutputValues(node, out);
          xOut += xLocal * nodeBasis;
        }
        ofs.write(reinterpret_cast<const char *>(&xOut), sizeof(double));
      }
    }
  }

  // DUMP Y
  {
    std::ofstream ofs(dataOutputDir / "y.fem", std::ios::binary);
    constexpr std::size_t direction = 1;
    for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
      for (std::size_t out = 0; out < kNumberOfOutputPointsPerCell; out++) {
        double yOut = 0.0;
        for (std::size_t node = 0; node < kNumberOfKinematicPointsPerCell;
             node++) {
          const std::size_t nodeIndex =
              getKinematicIndexFromCell(cell, node, direction);
          const double yLocal = x(nodeIndex);
          const double nodeBasis = kinematicBasisOutputValues(node, out);
          yOut += yLocal * nodeBasis;
        }
        ofs.write(reinterpret_cast<const char *>(&yOut), sizeof(double));
      }
    }
  }

  // DUMP U
  {
    std::ofstream ofs(dataOutputDir / "u.fem", std::ios::binary);
    constexpr std::size_t direction = 0;
    for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
      for (std::size_t out = 0; out < kNumberOfOutputPointsPerCell; out++) {
        double uOut = 0.0;
        for (std::size_t node = 0; node < kNumberOfKinematicPointsPerCell;
             node++) {
          const std::size_t nodeIndex =
              getKinematicIndexFromCell(cell, node, direction);
          const double uLocal = u(nodeIndex);
          const double nodeBasis = kinematicBasisOutputValues(node, out);
          uOut += uLocal * nodeBasis;
        }
        ofs.write(reinterpret_cast<const char *>(&uOut), sizeof(double));
      }
    }
  }

  // DUMP V
  {
    std::ofstream ofs(dataOutputDir / "v.fem", std::ios::binary);
    constexpr std::size_t direction = 1;
    for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
      for (std::size_t out = 0; out < kNumberOfOutputPointsPerCell; out++) {
        double vOut = 0.0;
        for (std::size_t node = 0; node < kNumberOfKinematicPointsPerCell;
             node++) {
          const std::size_t nodeIndex =
              getKinematicIndexFromCell(cell, node, direction);
          const double vLocal = u(nodeIndex);
          const double nodeBasis = kinematicBasisOutputValues(node, out);
          vOut += vLocal * nodeBasis;
        }
        ofs.write(reinterpret_cast<const char *>(&vOut), sizeof(double));
      }
    }
  }

  // DUMP VOLUME FRACTIONS, RHO, E, P
  {
    for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
      const std::string filenameVolFrac =
          "vol_frac" + std::to_string(material) + ".fem";
      const std::string filenameRho = "rho" + std::to_string(material) + ".fem";
      const std::string filenameE = "e" + std::to_string(material) + ".fem";
      const std::string filenameP = "p" + std::to_string(material) + ".fem";

      std::ofstream ofsVolFrac(dataOutputDir / filenameVolFrac,
                               std::ios::binary);
      std::ofstream ofsRho(dataOutputDir / filenameRho, std::ios::binary);
      std::ofstream ofsE(dataOutputDir / filenameE, std::ios::binary);
      std::ofstream ofsP(dataOutputDir / filenameP, std::ios::binary);

      for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
        for (std::size_t out = 0; out < kNumberOfOutputPointsPerCell; out++) {
          Eigen::Matrix2d jacobian;
          Eigen::Matrix2d jacobianInitial;
          jacobian.setZero();
          jacobianInitial.setZero();
          for (std::size_t node = 0; node < kNumberOfKinematicPointsPerCell;
               node++) {
            for (std::size_t direction = 0; direction < kSolverDimention;
                 direction++) {
              const std::size_t nodeIndex =
                  getKinematicIndexFromCell(cell, node, direction);
              const double xLocal = x(nodeIndex);
              const double xInitialLocal = xInitial(nodeIndex);
              const double nodeBasisdx =
                  kinematicBasisdxOutputValues(node, out);
              const double nodeBasisdy =
                  kinematicBasisdyOutputValues(node, out);
              jacobian(0, direction) += xLocal * nodeBasisdx;
              jacobian(1, direction) += xLocal * nodeBasisdy;
              jacobianInitial(0, direction) += xInitialLocal * nodeBasisdx;
              jacobianInitial(1, direction) += xInitialLocal * nodeBasisdy;
            }
          }
          const double jacobianDet = std::abs(jacobian.determinant());
          const double jacobianInitialDet =
              std::abs(jacobianInitial.determinant());

          double volFracInitialOut = 0.0;
          double rhoInitialOut = 0.0;
          double volFracOut = 0.0;

          for (std::size_t quad = 0; quad < kNumberOfQuadraturePointsPerCell;
               quad++) {
            const std::size_t quadIndex =
                getQuadIndexFromCell(cell, quad, material);
            const double volFracInitialLocal = volFracInitial(quadIndex);
            const double rhoInitialLocal = rhoInitial(quadIndex);
            const double volFracLocal = volFrac(quadIndex);
            const double quadBasis = quadBasisOutputValues(quad, out);

            volFracInitialOut += volFracInitialLocal * quadBasis;
            rhoInitialOut += rhoInitialLocal * quadBasis;
            volFracOut += volFracLocal * quadBasis;
          }

          const double rhoOut =
              volFracInitialOut > 0.0
                  ? (volFracInitialOut * rhoInitialOut * jacobianInitialDet) /
                        (volFracOut * jacobianDet)
                  : 0.0;

          double eOut = 0.0;

          for (std::size_t node = 0; node < kNumberOfThermodynamicPointsPerCell;
               node++) {
            const std::size_t nodeIndex =
                getThermodynamicIndexFromCell(cell, node, material);

            const double eLocal = e(nodeIndex);
            const double nodeBasis = thermodynamicBasisOutputValues(node, out);
            eOut += eLocal * nodeBasis;
          }

          auto eos = problem.eoses[material];

          double pOut = eos->getp(rhoOut, eOut);
          // const double pMin =
          //     (-1.0 / 3.0) * problem.plastic_params->equivalent_stress;
          // if (pOut < pMin) {
          //   pOut = pMin;
          // }

          ofsVolFrac.write(reinterpret_cast<const char *>(&volFracOut),
                           sizeof(double));
          ofsRho.write(reinterpret_cast<const char *>(&rhoOut), sizeof(double));
          ofsE.write(reinterpret_cast<const char *>(&eOut), sizeof(double));
          ofsP.write(reinterpret_cast<const char *>(&pOut), sizeof(double));
        }
      }
    }
  }

  dumpPlasticityThreshold(dataOutputDir);
  dumpDeviatorStressTensor(dataOutputDir);
  dumpEquivalentPlasticStrain(dataOutputDir);
}

void PlasticElasticMethod::calcStep() { RK2Step(); }

double PlasticElasticMethod::initHminCoeff(
    const std::size_t numberOfKinematicPointsPerCellPerDimention) {
  assert(numberOfKinematicPointsPerCellPerDimention >= 2);

  const std::size_t kinematicOrder =
      numberOfKinematicPointsPerCellPerDimention - 1;

  const std::vector<double> nodes = getLobattoNodes(kinematicOrder);
  assert(nodes.size() == numberOfKinematicPointsPerCellPerDimention);

  return nodes[1];
}

std::vector<double> PlasticElasticMethod::initQuadratureWeights(
    const std::size_t numberOfQuadPointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const std::size_t numberOfQuadPoints =
      std::pow(numberOfQuadPointsPerCellPerDimention, solverDimention);

  const std::vector<double> quadWeights1D =
      getLegendreWeights(numberOfQuadPointsPerCellPerDimention - 1);
  assert(quadWeights1D.size() == numberOfQuadPointsPerCellPerDimention);

  std::vector<double> output;
  output.reserve(numberOfQuadPoints);
  for (std::size_t i = 0; i < quadWeights1D.size(); i++) {
    for (std::size_t j = 0; j < quadWeights1D.size(); j++) {
      output.push_back(quadWeights1D[i] * quadWeights1D[j]);
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::initKinematicBasisQuadValues(
    const std::size_t numberOfKinematicPointsPerCellPerDimention,
    const std::size_t numberOfQuadraturePointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const std::vector<double> quadNodes1D =
      getLegendreNodes(numberOfQuadraturePointsPerCellPerDimention - 1);
  assert(quadNodes1D.size() == numberOfQuadraturePointsPerCellPerDimention);

  const std::size_t numberOfKinematicPointsPerCell =
      std::pow(numberOfKinematicPointsPerCellPerDimention, solverDimention);
  const std::size_t numberOfQuadraturePointsPerCell =
      std::pow(numberOfQuadraturePointsPerCellPerDimention, solverDimention);

  Eigen::MatrixXd output(numberOfKinematicPointsPerCell,
                         numberOfQuadraturePointsPerCell);
  output.setZero();

  const std::size_t lobattoOrder =
      numberOfKinematicPointsPerCellPerDimention - 1;
  for (std::size_t quad = 0; quad < numberOfQuadraturePointsPerCell; quad++) {
    const std::size_t quadi =
        quad / numberOfQuadraturePointsPerCellPerDimention;
    const std::size_t quadj =
        quad % numberOfQuadraturePointsPerCellPerDimention;

    const double x = quadNodes1D[quadi];
    const double y = quadNodes1D[quadj];

    for (std::size_t basis = 0; basis < numberOfKinematicPointsPerCell;
         basis++) {
      const std::size_t basisi =
          basis / numberOfKinematicPointsPerCellPerDimention;
      const std::size_t basisj =
          basis % numberOfKinematicPointsPerCellPerDimention;

      const double basisx = lobattoBasis1D(x, lobattoOrder, basisi);
      const double basisy = lobattoBasis1D(y, lobattoOrder, basisj);

      output(basis, quad) = basisx * basisy;
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::initKinematicBasisdxQuadValues(
    const std::size_t numberOfKinematicPointsPerCellPerDimention,
    const std::size_t numberOfQuadraturePointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const std::vector<double> quadNodes1D =
      getLegendreNodes(numberOfQuadraturePointsPerCellPerDimention - 1);
  assert(quadNodes1D.size() == numberOfQuadraturePointsPerCellPerDimention);

  const std::size_t numberOfKinematicPointsPerCell =
      std::pow(numberOfKinematicPointsPerCellPerDimention, solverDimention);
  const std::size_t numberOfQuadraturePointsPerCell =
      std::pow(numberOfQuadraturePointsPerCellPerDimention, solverDimention);

  Eigen::MatrixXd output(numberOfKinematicPointsPerCell,
                         numberOfQuadraturePointsPerCell);
  output.setZero();

  const std::size_t lobattoOrder =
      numberOfKinematicPointsPerCellPerDimention - 1;

  for (std::size_t quad = 0; quad < numberOfQuadraturePointsPerCell; quad++) {
    const std::size_t quadi =
        quad / numberOfQuadraturePointsPerCellPerDimention;
    const std::size_t quadj =
        quad % numberOfQuadraturePointsPerCellPerDimention;

    const double x = quadNodes1D[quadi];
    const double y = quadNodes1D[quadj];

    for (std::size_t basis = 0; basis < numberOfKinematicPointsPerCell;
         basis++) {
      const std::size_t basisi =
          basis / numberOfKinematicPointsPerCellPerDimention;
      const std::size_t basisj =
          basis % numberOfKinematicPointsPerCellPerDimention;

      const double basisx = lobattoBasis1Ddx(x, lobattoOrder, basisi);
      const double basisy = lobattoBasis1D(y, lobattoOrder, basisj);

      output(basis, quad) = basisx * basisy;
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::initKinematicBasisdyQuadValues(
    const std::size_t numberOfKinematicPointsPerCellPerDimention,
    const std::size_t numberOfQuadraturePointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const std::vector<double> quadNodes1D =
      getLegendreNodes(numberOfQuadraturePointsPerCellPerDimention - 1);
  assert(quadNodes1D.size() == numberOfQuadraturePointsPerCellPerDimention);

  const std::size_t numberOfKinematicPointsPerCell =
      std::pow(numberOfKinematicPointsPerCellPerDimention, solverDimention);
  const std::size_t numberOfQuadraturePointsPerCell =
      std::pow(numberOfQuadraturePointsPerCellPerDimention, solverDimention);

  Eigen::MatrixXd output(numberOfKinematicPointsPerCell,
                         numberOfQuadraturePointsPerCell);
  output.setZero();

  const std::size_t lobattoOrder =
      numberOfKinematicPointsPerCellPerDimention - 1;

  for (std::size_t quad = 0; quad < numberOfQuadraturePointsPerCell; quad++) {
    const std::size_t quadi =
        quad / numberOfQuadraturePointsPerCellPerDimention;
    const std::size_t quadj =
        quad % numberOfQuadraturePointsPerCellPerDimention;

    const double x = quadNodes1D[quadi];
    const double y = quadNodes1D[quadj];

    for (std::size_t basis = 0; basis < numberOfKinematicPointsPerCell;
         basis++) {
      const std::size_t basisi =
          basis / numberOfKinematicPointsPerCellPerDimention;
      const std::size_t basisj =
          basis % numberOfKinematicPointsPerCellPerDimention;

      const double basisx = lobattoBasis1D(x, lobattoOrder, basisi);
      const double basisy = lobattoBasis1Ddx(y, lobattoOrder, basisj);

      output(basis, quad) = basisx * basisy;
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::initThermodynamicBasisQuadValues(
    const std::size_t numberOfThermodynamicPointsPerCellPerDimention,
    const std::size_t numberOfQuadraturePointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const std::vector<double> quadNodes1D =
      getLegendreNodes(numberOfQuadraturePointsPerCellPerDimention - 1);

  const std::size_t numberOfThermodynamicPointsPerCell =
      std::pow(numberOfThermodynamicPointsPerCellPerDimention, solverDimention);
  const std::size_t numberOfQuadraturePointsPerCell =
      std::pow(numberOfQuadraturePointsPerCellPerDimention, solverDimention);

  Eigen::MatrixXd output(numberOfThermodynamicPointsPerCell,
                         numberOfQuadraturePointsPerCell);
  output.setZero();

  const std::size_t legendreOrder =
      numberOfThermodynamicPointsPerCellPerDimention - 1;
  for (std::size_t quad = 0; quad < numberOfQuadraturePointsPerCell; quad++) {
    const std::size_t quadi =
        quad / numberOfQuadraturePointsPerCellPerDimention;
    const std::size_t quadj =
        quad % numberOfQuadraturePointsPerCellPerDimention;

    const double x = quadNodes1D[quadi];
    const double y = quadNodes1D[quadj];

    for (std::size_t basis = 0; basis < numberOfThermodynamicPointsPerCell;
         basis++) {
      const std::size_t basisi =
          basis / numberOfThermodynamicPointsPerCellPerDimention;
      const std::size_t basisj =
          basis % numberOfThermodynamicPointsPerCellPerDimention;

      const double basisx = legendreBasis1D(x, legendreOrder, basisi);
      const double basisy = legendreBasis1D(y, legendreOrder, basisj);

      output(basis, quad) = basisx * basisy;
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::initKinematicBasisOutputValues(
    const std::size_t numberOfKinematicPointsPerCellPerDimention,
    const std::size_t numberOfOutputPointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const double di = 1.0 / numberOfOutputPointsPerCellPerDimention;
  std::vector<double> outputNodes1D(numberOfOutputPointsPerCellPerDimention);
  for (std::size_t i = 0; i < numberOfOutputPointsPerCellPerDimention; i++) {
    outputNodes1D[i] = 0.5 * di + i * di;
  }

  const std::size_t numberOfKinematicPointsPerCell =
      std::pow(numberOfKinematicPointsPerCellPerDimention, solverDimention);
  const std::size_t numberOfOutputPointsPerCell =
      std::pow(numberOfOutputPointsPerCellPerDimention, solverDimention);

  Eigen::MatrixXd output(numberOfKinematicPointsPerCell,
                         numberOfOutputPointsPerCell);
  output.setZero();

  const std::size_t lobattoOrder =
      numberOfKinematicPointsPerCellPerDimention - 1;
  for (std::size_t out = 0; out < numberOfOutputPointsPerCell; out++) {
    const std::size_t outi = out / numberOfOutputPointsPerCellPerDimention;
    const std::size_t outj = out % numberOfOutputPointsPerCellPerDimention;

    const double x = outputNodes1D[outi];
    const double y = outputNodes1D[outj];

    for (std::size_t node = 0; node < numberOfKinematicPointsPerCell; node++) {
      const std::size_t nodei =
          node / numberOfKinematicPointsPerCellPerDimention;
      const std::size_t nodej =
          node % numberOfKinematicPointsPerCellPerDimention;

      const double basisi = lobattoBasis1D(x, lobattoOrder, nodei);
      const double basisj = lobattoBasis1D(y, lobattoOrder, nodej);

      output(node, out) = basisi * basisj;
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::initKinematicBasisdxOutputValues(
    const std::size_t numberOfKinematicPointsPerCellPerDimention,
    const std::size_t numberOfOutputPointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const double di = 1.0 / numberOfOutputPointsPerCellPerDimention;
  std::vector<double> outputNodes1D(numberOfOutputPointsPerCellPerDimention);
  for (std::size_t i = 0; i < numberOfOutputPointsPerCellPerDimention; i++) {
    outputNodes1D[i] = 0.5 * di + i * di;
  }

  const std::size_t numberOfKinematicPointsPerCell =
      std::pow(numberOfKinematicPointsPerCellPerDimention, solverDimention);
  const std::size_t numberOfOutputPointsPerCell =
      std::pow(numberOfOutputPointsPerCellPerDimention, solverDimention);

  Eigen::MatrixXd output(numberOfKinematicPointsPerCell,
                         numberOfOutputPointsPerCell);
  output.setZero();

  const std::size_t lobattoOrder =
      numberOfKinematicPointsPerCellPerDimention - 1;
  for (std::size_t out = 0; out < numberOfOutputPointsPerCell; out++) {
    const std::size_t outi = out / numberOfOutputPointsPerCellPerDimention;
    const std::size_t outj = out % numberOfOutputPointsPerCellPerDimention;

    const double x = outputNodes1D[outi];
    const double y = outputNodes1D[outj];

    for (std::size_t node = 0; node < numberOfKinematicPointsPerCell; node++) {
      const std::size_t nodei =
          node / numberOfKinematicPointsPerCellPerDimention;
      const std::size_t nodej =
          node % numberOfKinematicPointsPerCellPerDimention;

      const double basisi = lobattoBasis1Ddx(x, lobattoOrder, nodei);
      const double basisj = lobattoBasis1D(y, lobattoOrder, nodej);

      output(node, out) = basisi * basisj;
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::initKinematicBasisdyOutputValues(
    const std::size_t numberOfKinematicPointsPerCellPerDimention,
    const std::size_t numberOfOutputPointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const double di = 1.0 / numberOfOutputPointsPerCellPerDimention;
  std::vector<double> outputNodes1D(numberOfOutputPointsPerCellPerDimention);
  for (std::size_t i = 0; i < numberOfOutputPointsPerCellPerDimention; i++) {
    outputNodes1D[i] = 0.5 * di + i * di;
  }

  const std::size_t numberOfKinematicPointsPerCell =
      std::pow(numberOfKinematicPointsPerCellPerDimention, solverDimention);
  const std::size_t numberOfOutputPointsPerCell =
      std::pow(numberOfOutputPointsPerCellPerDimention, solverDimention);

  Eigen::MatrixXd output(numberOfKinematicPointsPerCell,
                         numberOfOutputPointsPerCell);
  output.setZero();

  const std::size_t lobattoOrder =
      numberOfKinematicPointsPerCellPerDimention - 1;
  for (std::size_t out = 0; out < numberOfOutputPointsPerCell; out++) {
    const std::size_t outi = out / numberOfOutputPointsPerCellPerDimention;
    const std::size_t outj = out % numberOfOutputPointsPerCellPerDimention;

    const double x = outputNodes1D[outi];
    const double y = outputNodes1D[outj];

    for (std::size_t node = 0; node < numberOfKinematicPointsPerCell; node++) {
      const std::size_t nodei =
          node / numberOfKinematicPointsPerCellPerDimention;
      const std::size_t nodej =
          node % numberOfKinematicPointsPerCellPerDimention;

      const double basisi = lobattoBasis1D(x, lobattoOrder, nodei);
      const double basisj = lobattoBasis1Ddx(y, lobattoOrder, nodej);

      output(node, out) = basisi * basisj;
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::initThermodynamicBasisOutputValues(
    const std::size_t numberOfThermodynamicPointsPerCellPerDimention,
    const std::size_t numberOfOutputPointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const double di = 1.0 / numberOfOutputPointsPerCellPerDimention;
  std::vector<double> outputNodes1D(numberOfOutputPointsPerCellPerDimention);
  for (std::size_t i = 0; i < numberOfOutputPointsPerCellPerDimention; i++) {
    outputNodes1D[i] = 0.5 * di + i * di;
  }

  const std::size_t numberOfThermodynamicPointsPerCell =
      std::pow(numberOfThermodynamicPointsPerCellPerDimention, solverDimention);
  const std::size_t numberOfOutputPointsPerCell =
      std::pow(numberOfOutputPointsPerCellPerDimention, solverDimention);

  Eigen::MatrixXd output(numberOfThermodynamicPointsPerCell,
                         numberOfOutputPointsPerCell);
  output.setZero();

  const std::size_t legendreOrder =
      numberOfThermodynamicPointsPerCellPerDimention - 1;
  for (std::size_t out = 0; out < numberOfOutputPointsPerCell; out++) {
    const std::size_t outi = out / numberOfOutputPointsPerCellPerDimention;
    const std::size_t outj = out % numberOfOutputPointsPerCellPerDimention;

    const double x = outputNodes1D[outi];
    const double y = outputNodes1D[outj];

    for (std::size_t node = 0; node < numberOfThermodynamicPointsPerCell;
         node++) {
      const std::size_t nodei =
          node / numberOfThermodynamicPointsPerCellPerDimention;
      const std::size_t nodej =
          node % numberOfThermodynamicPointsPerCellPerDimention;

      const double basisi = legendreBasis1D(x, legendreOrder, nodei);
      const double basisj = legendreBasis1D(y, legendreOrder, nodej);

      output(node, out) = basisi * basisj;
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::initQuadBasisOutputValues(
    const std::size_t numberOfQuadraturePointsPerCellPerDimention,
    const std::size_t numberOfOutputPointsPerCellPerDimention,
    const std::size_t solverDimention) {
  assert(solverDimention == 2);

  const double di = 1.0 / numberOfOutputPointsPerCellPerDimention;
  std::vector<double> outputNodes1D(numberOfOutputPointsPerCellPerDimention);
  for (std::size_t i = 0; i < numberOfOutputPointsPerCellPerDimention; i++) {
    outputNodes1D[i] = 0.5 * di + i * di;
  }

  const std::size_t numberOfQuadraturePointsPerCell =
      std::pow(numberOfQuadraturePointsPerCellPerDimention, solverDimention);
  const std::size_t numberOfOutputPointsPerCell =
      std::pow(numberOfOutputPointsPerCellPerDimention, solverDimention);

  Eigen::MatrixXd output(numberOfQuadraturePointsPerCell,
                         numberOfOutputPointsPerCell);
  output.setZero();

  const std::size_t legendreOrder =
      numberOfQuadraturePointsPerCellPerDimention - 1;
  for (std::size_t out = 0; out < numberOfOutputPointsPerCell; out++) {
    const std::size_t outi = out / numberOfOutputPointsPerCellPerDimention;
    const std::size_t outj = out % numberOfOutputPointsPerCellPerDimention;

    const double x = outputNodes1D[outi];
    const double y = outputNodes1D[outj];

    for (std::size_t quad = 0; quad < numberOfQuadraturePointsPerCell; quad++) {
      const std::size_t quadi =
          quad / numberOfQuadraturePointsPerCellPerDimention;
      const std::size_t quadj =
          quad % numberOfQuadraturePointsPerCellPerDimention;

      const double basisi = legendreBasis1D(x, legendreOrder, quadi);
      const double basisj = legendreBasis1D(y, legendreOrder, quadj);

      output(quad, out) = basisi * basisj;
    }
  }

  return output;
}

void PlasticElasticMethod::initKinematicVectors() {
  assert(kNumberOfKinematicPointsPerCellPerDimention >= 1);

  const double celldx = (problem.xmax - problem.xmin) / xCells;
  const double celldy = (problem.ymax - problem.ymin) / yCells;

  switch (problem.dimension) {
    case ProblemDimension::e1D:
      l0 = celldx / kOrder;
      break;
    case ProblemDimension::e2D:
      l0 = std::sqrt(celldx * celldy) / kOrder;
      break;
  }

  const std::size_t kinematicOrder =
      kNumberOfKinematicPointsPerCellPerDimention - 1;

  const std::vector<double> kinematicNodes1D = getLobattoNodes(kinematicOrder);
  assert(kinematicNodes1D.size() ==
         kNumberOfKinematicPointsPerCellPerDimention);

  x.resize(kSolverDimention * kNumberOfKinematicPointsTotal);
  u.resize(kSolverDimention * kNumberOfKinematicPointsTotal);

  const std::size_t imax = xCells * kinematicOrder + 1;
  const std::size_t jmax = yCells * kinematicOrder + 1;
  for (std::size_t i = 0; i < imax; i++) {
    for (std::size_t j = 0; j < jmax; j++) {
      const double localx = kinematicNodes1D[i % kinematicOrder];
      const double localy = kinematicNodes1D[j % kinematicOrder];

      const std::size_t celli = i / kinematicOrder;
      const std::size_t cellj = j / kinematicOrder;
      const double xij = problem.xmin + celldx * (celli + localx);
      const double yij = problem.ymin + celldy * (cellj + localy);

      x(kSolverDimention * (i * jmax + j)) = xij;
      x(kSolverDimention * (i * jmax + j) + 1) = yij;

      u(kSolverDimention * (i * jmax + j)) = problem.uInitializer(xij, yij);
      u(kSolverDimention * (i * jmax + j) + 1) = problem.vInitializer(xij, yij);
    }
  }

  xInitial = x;
}

void PlasticElasticMethod::initQuadVectors() {
  initVolFracVector();
  initRhoVector();
}

void PlasticElasticMethod::initVolFracVector() {
  assert(kNumberOfQuadraturePointsPerCellPerDimention >= 1);

  const double celldx = (problem.xmax - problem.xmin) / xCells;
  const double celldy = (problem.ymax - problem.ymin) / yCells;

  const std::size_t quadratureOrder =
      kNumberOfQuadraturePointsPerCellPerDimention - 1;

  const std::vector<double> quadratureNodes1D =
      getLegendreNodes(quadratureOrder);
  assert(quadratureNodes1D.size() ==
         kNumberOfQuadraturePointsPerCellPerDimention);

  volFrac.resize(kNumberOfMaterials * kNumberOfQuadraturePointsTotal);
  volFrac.setZero();

  assert(kOrder >= 1);
  const std::size_t volFracInterpolOrder = kOrder - 1;
  const double di = volFracInterpolOrder > 0 ? 1.0 / volFracInterpolOrder : 0.0;
  for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
    const std::size_t celli = cell / yCells;
    const std::size_t cellj = cell % yCells;

    for (std::size_t i = 0; i <= volFracInterpolOrder; i++) {
      const double localx = i * di;
      const double xij = problem.xmin + celldx * (celli + localx);
      for (std::size_t j = 0; j <= volFracInterpolOrder; j++) {
        const double localy = j * di;

        const double yij = problem.ymin + celldy * (cellj + localy);

        const std::vector<double> volFracs =
            problem.volFractionInitializer(xij, yij);
        assert(volFracs.size() == kNumberOfMaterials);

        for (std::size_t quad = 0; quad < kNumberOfQuadraturePointsPerCell;
             quad++) {
          const std::size_t quadi =
              quad / kNumberOfQuadraturePointsPerCellPerDimention;
          const std::size_t quadj =
              quad % kNumberOfQuadraturePointsPerCellPerDimention;

          const double quadx = quadratureNodes1D[quadi];
          const double quady = quadratureNodes1D[quadj];

          const double basisi =
              bernsteinBasis1D(quadx, volFracInterpolOrder, i);
          const double basisj =
              bernsteinBasis1D(quady, volFracInterpolOrder, j);

          for (std::size_t material = 0; material < kNumberOfMaterials;
               material++) {
            const std::size_t quadIndex =
                getQuadIndexFromCell(cell, quad, material);

            volFrac(quadIndex) += volFracs[material] * basisi * basisj;
          }
        }
      }
    }
  }

  volFracInitial = volFrac;
}

void PlasticElasticMethod::initRhoVector() {
  assert(kNumberOfQuadraturePointsPerCellPerDimention >= 1);

  const double celldx = (problem.xmax - problem.xmin) / xCells;
  const double celldy = (problem.ymax - problem.ymin) / yCells;

  const std::size_t quadratureOrder =
      kNumberOfQuadraturePointsPerCellPerDimention - 1;

  const std::vector<double> quadratureNodes1D =
      getLegendreNodes(quadratureOrder);
  assert(quadratureNodes1D.size() ==
         kNumberOfQuadraturePointsPerCellPerDimention);

  rhoInitial.resize(kNumberOfMaterials * kNumberOfQuadraturePointsTotal);
  rhoInitial.setZero();

  for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
    const std::size_t celli = cell / yCells;
    const std::size_t cellj = cell % yCells;

    for (std::size_t quad = 0; quad < kNumberOfQuadraturePointsPerCell;
         quad++) {
      // const std::size_t quadi =
      //     quad / kNumberOfQuadraturePointsPerCellPerDimention;
      // const std::size_t quadj =
      //     quad % kNumberOfQuadraturePointsPerCellPerDimention;

      // const double localx = quadratureNodes1D[quadi];
      // const double localy = quadratureNodes1D[quadj];
      const double localx = 0.5;
      const double localy = 0.5;

      const double xij = problem.xmin + celldx * (celli + localx);
      const double yij = problem.ymin + celldy * (cellj + localy);

      const std::vector<double> rhos = problem.rhoInitializer(xij, yij);
      assert(rhos.size() == kNumberOfMaterials);

      for (std::size_t material = 0; material < kNumberOfMaterials;
           material++) {
        const std::size_t quadIndex =
            getQuadIndexFromCell(cell, quad, material);

        if (volFrac(quadIndex) == 0.0) {
          continue;
        }

        rhoInitial(quadIndex) = rhos[material];
      }
    }
  }
}

void PlasticElasticMethod::initThermodynamicVector() {
  assert(kNumberOfThermodynamicPointsPerCellPerDimention >= 1);

  const double celldx = (problem.xmax - problem.xmin) / xCells;
  const double celldy = (problem.ymax - problem.ymin) / yCells;

  const std::size_t thermodynamicOrder =
      kNumberOfThermodynamicPointsPerCellPerDimention - 1;

  const std::vector<double> thermoNodes1D =
      getLegendreNodes(thermodynamicOrder);
  assert(thermoNodes1D.size() ==
         kNumberOfThermodynamicPointsPerCellPerDimention);

  e.resize(kNumberOfMaterials * kNumberOfThermodynamicPointsTotal);
  e.setZero();

  for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
    const std::size_t celli = cell / yCells;
    const std::size_t cellj = cell % yCells;

    for (std::size_t node = 0; node < kNumberOfThermodynamicPointsPerCell;
         node++) {
      // const double localx = thermoNodes1D[nodei];
      // const double localy = thermoNodes1D[nodej];
      const double localx = 0.5;
      const double localy = 0.5;

      const double xij = problem.xmin + celldx * (celli + localx);
      const double yij = problem.ymin + celldy * (cellj + localy);

      const std::vector<double> rhos = problem.rhoInitializer(xij, yij);
      const std::vector<double> ps = problem.pInitializer(xij, yij);

      for (std::size_t material = 0; material < kNumberOfMaterials;
           material++) {
        bool bVolFracIsZero = true;
        for (std::size_t quad = 0; quad < kNumberOfQuadraturePointsPerCell;
             quad++) {
          const std::size_t quadIndex =
              getQuadIndexFromCell(cell, quad, material);
          if (volFrac(quadIndex) != 0.0) {
            bVolFracIsZero = false;
            break;
          }
        }
        if (bVolFracIsZero) {
          continue;
        }
        const double rhoLocal = rhos[material];
        const double pLocal = ps[material];
        const double eLocal = problem.eoses[material]->gete(rhoLocal, pLocal);
        const std::size_t thermoIndex =
            getThermodynamicIndexFromCell(cell, node, material);

        e(thermoIndex) = eLocal;
      }
    }
  }
}

void PlasticElasticMethod::initKinematicMassMatrix() {
  kinematicMassMatrix.resize(kSolverDimention * kNumberOfKinematicPointsTotal,
                             kSolverDimention * kNumberOfKinematicPointsTotal);

  calcKinematicMassMatrix();
}

void PlasticElasticMethod::initThermoMassMatrixInv() {
  thermoMassMatrixInv.resize(
      kNumberOfMaterials * kNumberOfThermodynamicPointsTotal,
      kNumberOfMaterials * kNumberOfThermodynamicPointsTotal);

  thermoMassMatrixInv.reserve(Eigen::VectorXi::Constant(
      kNumberOfMaterials * kNumberOfThermodynamicPointsTotal,
      kNumberOfThermodynamicPointsPerCell));

  calcThermoMassMatrixInv();

  thermoMassMatrixInv.makeCompressed();
}

void PlasticElasticMethod::initForceMatrix() {
  forceMatrix.resize(kSolverDimention * kNumberOfKinematicPointsTotal,
                     kNumberOfMaterials * kNumberOfThermodynamicPointsTotal);
  forceMatrix.reserve(Eigen::VectorXi::Constant(
      kNumberOfMaterials * kNumberOfThermodynamicPointsTotal,
      kSolverDimention * kNumberOfKinematicPointsPerCell));
}

void PlasticElasticMethod::initVolFracRateVector() {
  volFracRate.resize(kNumberOfMaterials * kNumberOfQuadraturePointsTotal);
  volFracRate.setZero();
}

void PlasticElasticMethod::initEnergyExchangeRateVector() {
  energyExchangeRate.resize(kNumberOfMaterials *
                            kNumberOfThermodynamicPointsTotal);
  volFracRate.setZero();
}

void PlasticElasticMethod::initDeviatorStressTensorVectors() {
  const std::size_t N = kNumberOfMaterials * kNumberOfThermodynamicPointsTotal;

  stress_deviator_00.resize(N);
  stress_deviator_01.resize(N);
  stress_deviator_10.resize(N);
  stress_deviator_11.resize(N);

  stress_deviator_00.setZero();
  stress_deviator_01.setZero();
  stress_deviator_10.setZero();
  stress_deviator_11.setZero();

  // const std::size_t i_max = N / 2;
  // for (std::size_t i = 0; i < i_max; i++) {
  //   stress_deviator_00(i) = -0.2e9;
  //   stress_deviator_11(i) = 0.1e9;
  // }
}

void PlasticElasticMethod::initDeviatorStressRateVectors() {
  const std::size_t N = kNumberOfMaterials * kNumberOfThermodynamicPointsTotal;

  deviator_rate_00.resize(N);
  deviator_rate_01.resize(N);
  deviator_rate_10.resize(N);
  deviator_rate_11.resize(N);

  deviator_rate_00.setZero();
  deviator_rate_01.setZero();
  deviator_rate_10.setZero();
  deviator_rate_11.setZero();
}

void PlasticElasticMethod::initEquivalentPlasticStrainVector() {
  const std::size_t N = kNumberOfMaterials * kNumberOfThermodynamicPointsTotal;

  equivalent_plastic_strain.resize(N);
  equivalent_plastic_strain.setZero();
}

void PlasticElasticMethod::initKinematicSolver() {
  kinematicMassMatrixSolver.compute(kinematicMassMatrix);
  if (kinematicMassMatrixSolver.info() != Eigen::Success) {
    throw std::runtime_error("Failed to init kinematic solver");
  }
}

void PlasticElasticMethod::calcKinematicMassMatrix() {
  const std::size_t kinematicOrder =
      kNumberOfKinematicPointsPerCellPerDimention - 1;

  kinematicMassMatrix.setZero();
  kinematicMassMatrix.reserve(Eigen::VectorXi::Constant(
      kSolverDimention * kNumberOfKinematicPointsTotal,
      kSolverDimention * (2 * kinematicOrder + 1) * (2 * kinematicOrder + 1)));

#pragma omp parallel for
  for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
    const Eigen::MatrixXd cellKinematicMass = quadKinematicCellMass(cell);
    assert(cellKinematicMass.rows() == static_cast<Eigen::MatrixXd::Index>(
                                           kNumberOfKinematicPointsPerCell));
    assert(cellKinematicMass.cols() == static_cast<Eigen::MatrixXd::Index>(
                                           kNumberOfKinematicPointsPerCell));

    for (std::size_t nodei = 0; nodei < kNumberOfKinematicPointsPerCell;
         nodei++) {
      for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
           nodej++) {
        for (std::size_t direction = 0; direction < kSolverDimention;
             direction++) {
          const std::size_t nodeiIndex =
              getKinematicIndexFromCell(cell, nodei, direction);
          const std::size_t nodejIndex =
              getKinematicIndexFromCell(cell, nodej, direction);

          // #pragma omp critical(kinematicMassMatrixUpdate)
          {
            kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) +=
                cellKinematicMass(nodei, nodej);
          }
        }
      }
    }
  }
  resolveBoundaryKinematicMassMatrix();
  kinematicMassMatrix.makeCompressed();
}

void PlasticElasticMethod::calcThermoMassMatrixInv() {
#pragma omp parallel for
  for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
    for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
      const Eigen::MatrixXd cellThermoMass = quadThermoCellMass(cell, material);

      const double cellThermoMassDet = cellThermoMass.determinant();

      Eigen::MatrixXd cellThermoMassInv(kNumberOfThermodynamicPointsPerCell,
                                        kNumberOfThermodynamicPointsPerCell);
      cellThermoMassInv.setZero();

      if (cellThermoMassDet != 0) {
        cellThermoMassInv = cellThermoMass.inverse();
      }

      for (std::size_t nodei = 0; nodei < kNumberOfThermodynamicPointsPerCell;
           nodei++) {
        for (std::size_t nodej = 0; nodej < kNumberOfThermodynamicPointsPerCell;
             nodej++) {
          const std::size_t nodeiIndex =
              getThermodynamicIndexFromCell(cell, nodei, material);
          const std::size_t nodejIndex =
              getThermodynamicIndexFromCell(cell, nodej, material);

          // #pragma omp critical(thermoMassMatrixInvUpdate)
          {
            thermoMassMatrixInv.coeffRef(nodeiIndex, nodejIndex) =
                cellThermoMassInv(nodei, nodej);
          }
        }
      }
    }
  }
}

void PlasticElasticMethod::calcForceMatrix(
    const Eigen::VectorXd &xCalc, const Eigen::VectorXd &uCalc,
    const Eigen::VectorXd &volFracCalc, const Eigen::VectorXd &eCalc,
    const Eigen::VectorXd &deviatorStressTensor_00,
    const Eigen::VectorXd &deviatorStressTensor_01,
    const Eigen::VectorXd &deviatorStressTensor_10,
    const Eigen::VectorXd &deviatorStressTensor_11) {
  assert(xCalc.size() == x.size());
  assert(uCalc.size() == u.size());
  assert(volFracCalc.size() == volFrac.size());
  assert(eCalc.size() == e.size());

  energyExchangeRate.setZero();
  deviator_rate_00.setZero();
  deviator_rate_01.setZero();
  deviator_rate_10.setZero();
  deviator_rate_11.setZero();

#pragma omp parallel for
  for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
    const Eigen::MatrixXd cellForce =
        quadForceMatrix(cell, xCalc, uCalc, volFracCalc, eCalc,
                        deviatorStressTensor_00, deviatorStressTensor_01,
                        deviatorStressTensor_10, deviatorStressTensor_11);
    for (std::size_t nodeKinematic = 0;
         nodeKinematic < kNumberOfKinematicPointsPerCell; nodeKinematic++) {
      for (std::size_t nodeThermo = 0;
           nodeThermo < kNumberOfThermodynamicPointsPerCell; nodeThermo++) {
        for (std::size_t direction = 0; direction < kSolverDimention;
             direction++) {
          for (std::size_t material = 0; material < kNumberOfMaterials;
               material++) {
            const std::size_t kinematicIndexLocal =
                getKinematicIndexLocal(nodeKinematic, direction);
            const std::size_t thermoIndexLocal =
                getThermoIndexLocal(nodeThermo, material);

            const std::size_t kinematicIndex =
                getKinematicIndexFromCell(cell, nodeKinematic, direction);
            const std::size_t thermoIndex =
                getThermodynamicIndexFromCell(cell, nodeThermo, material);

            // #pragma omp critical(forceMatrixUpdate)
            {
              forceMatrix.coeffRef(kinematicIndex, thermoIndex) =
                  cellForce(kinematicIndexLocal, thermoIndexLocal);
            }
          }
        }
      }
    }
  }
  forceMatrix.makeCompressed();
}

void PlasticElasticMethod::radialReturn(
    Eigen::VectorXd &deviatorStressTensor_00_calc,
    Eigen::VectorXd &deviatorStressTensor_01_calc,
    Eigen::VectorXd &deviatorStressTensor_10_calc,
    Eigen::VectorXd &deviatorStressTensor_11_calc,
    Eigen::VectorXd &equivalentPlasticStrain_calc) {
#pragma omp parallel for
  for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
    for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
      for (std::size_t thermo_node = 0;
           thermo_node < kNumberOfThermodynamicPointsPerCell; thermo_node++) {
        const std::size_t thermo_index =
            getThermodynamicIndexFromCell(cell, thermo_node, material);

        const double s00_local = deviatorStressTensor_00_calc(thermo_index);
        const double s01_local = deviatorStressTensor_01_calc(thermo_index);
        const double s10_local = deviatorStressTensor_10_calc(thermo_index);
        const double s11_local = deviatorStressTensor_11_calc(thermo_index);

        const double s22_local = -s00_local - s11_local;
        const double stress_measure =
            std::sqrt(1.5 * (s00_local * s00_local + s01_local * s01_local +
                             s10_local * s10_local + s11_local * s11_local +
                             s22_local * s22_local));

        const double equivalent_stress =
            problem.plastic_params->equivalent_stress;

        if (stress_measure <= equivalent_stress) {
          continue;
        }
        const double shear_modulus = problem.plastic_params->shear_modulus;

        const double equiv_plastic_strain_local =
            equivalentPlasticStrain_calc(thermo_index);

        const double strain_coeff = equivalent_stress / stress_measure;

        const double s00_local_new = s00_local * strain_coeff;
        const double s01_local_new = s01_local * strain_coeff;
        const double s10_local_new = s10_local * strain_coeff;
        const double s11_local_new = s11_local * strain_coeff;

        const double equiv_plastic_strain_local_new =
            equiv_plastic_strain_local +
            (stress_measure - equivalent_stress) / (3.0 * shear_modulus);

        deviatorStressTensor_00_calc(thermo_index) = s00_local_new;
        deviatorStressTensor_01_calc(thermo_index) = s01_local_new;
        deviatorStressTensor_10_calc(thermo_index) = s10_local_new;
        deviatorStressTensor_11_calc(thermo_index) = s11_local_new;

        equivalentPlasticStrain_calc(thermo_index) =
            equiv_plastic_strain_local_new;
      }
    }
  }
}

Eigen::MatrixXd PlasticElasticMethod::quadKinematicCellMass(
    const std::size_t cell) {
  Eigen::MatrixXd output(kNumberOfKinematicPointsPerCell,
                         kNumberOfKinematicPointsPerCell);
  output.setZero();

  for (std::size_t quad = 0; quad < kNumberOfQuadraturePointsPerCell; quad++) {
    const double quadWeight = quadWeights[quad];

    const Eigen::Matrix2d jacobian = calcJacobian(cell, quad, this->x);
    const double jacobianDet = std::abs(jacobian.determinant());

    double rhoLocal = 0.0;
    for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
      const std::size_t quadIndex = getQuadIndexFromCell(cell, quad, material);
      const double volFracMat = volFrac(quadIndex);
      const double rhoMat = rhoInitial(quadIndex);

      rhoLocal += volFracMat * rhoMat;
    }

    assert(rhoLocal > 0.0);

    for (std::size_t nodei = 0; nodei < kNumberOfKinematicPointsPerCell;
         nodei++) {
      const double basisi = kinematicBasisQuadValues(nodei, quad);
      for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
           nodej++) {
        const double basisj = kinematicBasisQuadValues(nodej, quad);

        output(nodei, nodej) +=
            quadWeight * rhoLocal * basisi * basisj * jacobianDet;
      }
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::quadThermoCellMass(
    const std::size_t cell, const std::size_t material) {
  Eigen::MatrixXd output(kNumberOfThermodynamicPointsPerCell,
                         kNumberOfThermodynamicPointsPerCell);
  output.setZero();

  for (std::size_t quad = 0; quad < kNumberOfQuadraturePointsPerCell; quad++) {
    const double quadWeight = quadWeights[quad];

    const Eigen::Matrix2d jacobian = calcJacobian(cell, quad, this->x);
    const double jacobianDet = std::abs(jacobian.determinant());

    const std::size_t quadIndex = getQuadIndexFromCell(cell, quad, material);
    const double volFracLocal = volFrac(quadIndex);
    const double rhoLocal = rhoInitial(quadIndex);

    for (std::size_t nodei = 0; nodei < kNumberOfThermodynamicPointsPerCell;
         nodei++) {
      const double basisi = thermodynamicBasisQuadValues(nodei, quad);
      for (std::size_t nodej = 0; nodej < kNumberOfThermodynamicPointsPerCell;
           nodej++) {
        const double basisj = thermodynamicBasisQuadValues(nodej, quad);

        output(nodei, nodej) += quadWeight * volFracLocal * rhoLocal * basisi *
                                basisj * jacobianDet;
      }
    }
  }

  return output;
}

Eigen::MatrixXd PlasticElasticMethod::quadForceMatrix(
    const std::size_t cell, const Eigen::VectorXd &xCalc,
    const Eigen::VectorXd &uCalc, const Eigen::VectorXd &volFracCalc,
    const Eigen::VectorXd &eCalc,
    const Eigen::VectorXd &deviatorStressTensor_00,
    const Eigen::VectorXd &deviatorStressTensor_01,
    const Eigen::VectorXd &deviatorStressTensor_10,
    const Eigen::VectorXd &deviatorStressTensor_11) {
  Eigen::MatrixXd output(
      kSolverDimention * kNumberOfKinematicPointsPerCell,
      kNumberOfMaterials * kNumberOfThermodynamicPointsPerCell);
  output.setZero();

  for (std::size_t quad = 0; quad < kNumberOfQuadraturePointsPerCell; quad++) {
    const double quadWeight = quadWeights[quad];

    const Eigen::Matrix2d jacobian = calcJacobian(cell, quad, xCalc);
    const Eigen::Matrix2d jacobianInv = jacobian.inverse();
    const double jacobianDet = std::abs(jacobian.determinant());

    Eigen::JacobiSVD<Eigen::Matrix2d> svd(jacobian);
    const double hmin = hminCoeff * svd.singularValues().minCoeff();

    std::vector<double> volFracs(kNumberOfMaterials, 0.0);
    std::vector<double> rhos(kNumberOfMaterials, 0.0);
    std::vector<double> ps(kNumberOfMaterials, 0.0);
    std::vector<double> soundSpeeds(kNumberOfMaterials, 0.0);
    double velocityScalarGrad = 0.0;

    for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
      const std::size_t quadIndex = getQuadIndexFromCell(cell, quad, material);
      const double volFracMat = volFracCalc(quadIndex);

      if (volFracMat == 0.0) {
        continue;
      }

      double soundSpeed = 0.0;
      double rhoLocal = 0.0;
      double pLocal = 0.0;
      double maxViscosityCoeff = 0.0;

      Eigen::Matrix2d stressTensor = calcStressTensor(
          cell, quad, material, soundSpeed, rhoLocal, pLocal, maxViscosityCoeff,
          velocityScalarGrad, jacobian, jacobianDet, jacobianInv, uCalc,
          volFracCalc, eCalc, deviatorStressTensor_00, deviatorStressTensor_01,
          deviatorStressTensor_10, deviatorStressTensor_11);

      for (std::size_t nodeThermo = 0;
           nodeThermo < kNumberOfThermodynamicPointsPerCell; nodeThermo++) {
        const std::size_t thermoIndex =
            getThermoIndexLocal(nodeThermo, material);
        const double basisThermo =
            thermodynamicBasisQuadValues(nodeThermo, quad);

        for (std::size_t nodeKinematic = 0;
             nodeKinematic < kNumberOfKinematicPointsPerCell; nodeKinematic++) {
          Eigen::Vector2d gradBasis{
              kinematicBasisdxQuadValues(nodeKinematic, quad),
              kinematicBasisdyQuadValues(nodeKinematic, quad)};
          gradBasis = jacobianInv * gradBasis;

          for (std::size_t direction = 0; direction < kSolverDimention;
               direction++) {
            const std::size_t kinematicIndex =
                getKinematicIndexLocal(nodeKinematic, direction);

            const double scalarProd =
                stressTensor(0, direction) * gradBasis(0) +
                stressTensor(1, direction) * gradBasis(1);

            output(kinematicIndex, thermoIndex) += quadWeight * scalarProd *
                                                   basisThermo * volFracMat *
                                                   jacobianDet;
          }
        }
      }

      volFracs[material] = volFracMat;
      rhos[material] = rhoLocal;
      ps[material] = pLocal;
      soundSpeeds[material] = soundSpeed;
      calcTau(hmin, soundSpeed, rhoLocal, maxViscosityCoeff);
    }

    const std::vector<double> rates = problem.matClosure->calcVolFracRates(
        volFracs, rhos, ps, soundSpeeds, cell, hmin, dt);
    assert(rates.size() == kNumberOfMaterials);

    for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
      const std::size_t quadIndex = getQuadIndexFromCell(cell, quad, material);
      volFracRate(quadIndex) = rates[material];
    }

    const double commonP = problem.matClosure->calcExchangePressure(
        volFracs, rhos, ps, soundSpeeds, rates, velocityScalarGrad, dt);

    for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
      const double volFracMat = volFracs[material];

      if (volFracMat == 0.0) {
        continue;
      }

      const double volFracMatRate = rates[material];

      for (std::size_t nodeThermo = 0;
           nodeThermo < kNumberOfThermodynamicPointsPerCell; nodeThermo++) {
        const std::size_t thermoIndex =
            getThermodynamicIndexFromCell(cell, nodeThermo, material);
        const double basisThermo =
            thermodynamicBasisQuadValues(nodeThermo, quad);

        const double value =
            quadWeight * commonP * volFracMatRate * basisThermo * jacobianDet;
        energyExchangeRate(thermoIndex) += value;
      }
    }
  }

  return output;
}

Eigen::Matrix2d PlasticElasticMethod::calcStressTensor(
    const std::size_t cell, const std::size_t quad, const std::size_t material,
    double &soundSpeed, double &rhoLocal, double &pLocal,
    double &maxViscosityCoeff, double &velocityScalarGrad,
    const Eigen::Matrix2d &jacobian, const double jacobianDet,
    const Eigen::Matrix2d &jacobianInv, const Eigen::VectorXd &uCalc,
    const Eigen::VectorXd &volFracCalc, const Eigen::VectorXd &eCalc,
    const Eigen::VectorXd &deviatorStressTensor_00_calc,
    const Eigen::VectorXd &deviatorStressTensor_01_calc,
    const Eigen::VectorXd &deviatorStressTensor_10_calc,
    const Eigen::VectorXd &deviatorStressTensor_11_calc) {
  Eigen::Matrix2d output;
  output.setZero();

  double eLocal = 0.0;
  for (std::size_t thermoNode = 0;
       thermoNode < kNumberOfThermodynamicPointsPerCell; thermoNode++) {
    const std::size_t thermoIndex =
        getThermodynamicIndexFromCell(cell, thermoNode, material);
    const double eNode = eCalc(thermoIndex);
    const double thermoBasis = thermodynamicBasisQuadValues(thermoNode, quad);
    eLocal += eNode * thermoBasis;
  }

  const Eigen::Matrix2d jacobianInitial = calcJacobian(cell, quad, xInitial);
  const double jacobianInitialDet = std::abs(jacobianInitial.determinant());
  assert(jacobianDet > 0.0);
  // const double jacobianDet = std::abs(jacobian.determinant());

  const std::size_t quadIndex = getQuadIndexFromCell(cell, quad, material);
  const double rhoInitialMat = rhoInitial(quadIndex);
  const double volFracInitialMat = volFracInitial(quadIndex);
  const double volFracMat = volFracCalc(quadIndex);

  rhoLocal = volFracInitialMat * rhoInitialMat * jacobianInitialDet /
             (jacobianDet * volFracMat);

  const auto eos = problem.eoses[material];
  pLocal = eos->getp(rhoLocal, eLocal);

  // const double pMin = (-1.0 / 3.0) *
  // problem.plastic_params->equivalent_stress; if (pLocal < pMin) {
  //   pLocal = pMin;
  // }
  soundSpeed = eos->getc(rhoLocal, pLocal);

  output = -pLocal * Eigen::Matrix2d::Identity();
  output += calcArtificialViscosity(
      cell, quad, soundSpeed, rhoLocal, maxViscosityCoeff, velocityScalarGrad,
      jacobian, jacobianInv, jacobianInitial, uCalc);
  output += calcDeviatorTensor(
      cell, quad, material, rhoLocal, jacobianDet, jacobianInv, uCalc,
      deviatorStressTensor_00_calc, deviatorStressTensor_01_calc,
      deviatorStressTensor_10_calc, deviatorStressTensor_11_calc);

  return output;
}

Eigen::Matrix2d PlasticElasticMethod::calcArtificialViscosity(
    const std::size_t cell, const std::size_t quad, const double soundSpeed,
    const double rhoLocal, double &maxViscosityCoeff,
    double &velocityScalarGrad, const Eigen::Matrix2d &jacobian,
    const Eigen::Matrix2d &jacobianInv, const Eigen::Matrix2d &jacobianInitial,
    const Eigen::VectorXd &uCalc) {
  Eigen::Matrix2d output;
  output.setZero();

  Eigen::Matrix2d velocityGrad;
  velocityGrad.setZero();

  for (std::size_t kinematicNode = 0;
       kinematicNode < kNumberOfKinematicPointsPerCell; kinematicNode++) {
    const double kinematicBasisdx =
        kinematicBasisdxQuadValues(kinematicNode, quad);
    const double kinematicBasisdy =
        kinematicBasisdyQuadValues(kinematicNode, quad);

    for (std::size_t direction = 0; direction < kSolverDimention; direction++) {
      const std::size_t kinematicIndex =
          getKinematicIndexFromCell(cell, kinematicNode, direction);
      const double uNode = uCalc(kinematicIndex);
      velocityGrad(0, direction) += uNode * kinematicBasisdx;
      velocityGrad(1, direction) += uNode * kinematicBasisdy;
    }
  }
  velocityGrad = jacobianInv * velocityGrad;
  velocityScalarGrad = velocityGrad(0, 0) + velocityGrad(1, 1);

  Eigen::Matrix2d symVelocityGrad =
      0.5 * (velocityGrad + velocityGrad.transpose());

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenSolver(symVelocityGrad);

  const double e1 = eigenSolver.eigenvalues()(0);
  const double e2 = eigenSolver.eigenvalues()(1);

  Eigen::Vector2d v1 = eigenSolver.eigenvectors().col(0);
  Eigen::Vector2d v2 = eigenSolver.eigenvectors().col(1);

  const double viscosityCoeff1 = calcViscosityCoeff(
      soundSpeed, rhoLocal, e1, v1, velocityGrad, jacobian, jacobianInitial);
  const double viscosityCoeff2 = calcViscosityCoeff(
      soundSpeed, rhoLocal, e2, v2, velocityGrad, jacobian, jacobianInitial);

  maxViscosityCoeff = std::max(viscosityCoeff1, viscosityCoeff2);

  const Eigen::Matrix2d v1Tensor{{v1(0) * v1(0), v1(0) * v1(1)},
                                 {v1(1) * v1(0), v1(1) * v1(1)}};
  const Eigen::Matrix2d v2Tensor{{v2(0) * v2(0), v2(0) * v2(1)},
                                 {v2(1) * v2(0), v2(1) * v2(1)}};

  output += viscosityCoeff1 * e1 * v1Tensor;
  output += viscosityCoeff2 * e2 * v2Tensor;

  return output;
}

double PlasticElasticMethod::calcViscosityCoeff(
    const double soundSpeed, const double rhoLocal, const double eigenvalue,
    const Eigen::Vector2d &eigenvector, const Eigen::Matrix2d &velocityGrad,
    const Eigen::Matrix2d &jacobian, const Eigen::Matrix2d &jacobianInitial) {
  const double velScalarGrad =
      std::abs(velocityGrad(0, 0) + velocityGrad(1, 1));
  const double velgradNorm = velocityGrad.norm();

  const double psi0 = velgradNorm != 0.0 ? velScalarGrad / velgradNorm : 0.0;

  const double psi1 = eigenvalue < 0 ? 1.0 : 0.0;

  const Eigen::Matrix2d jacobianInitialInv = jacobianInitial.inverse();
  const Eigen::Matrix2d jacobianMapping = jacobianInitialInv * jacobian;

  const double localLengthScale =
      l0 * (jacobianMapping * eigenvector).norm() / eigenvector.norm();

  const double quadraticCoeff =
      q2 * localLengthScale * localLengthScale * std::abs(eigenvalue);
  const double linearCoeff = q1 * psi0 * psi1 * localLengthScale * soundSpeed;
  const double output = rhoLocal * (quadraticCoeff + linearCoeff);
  return output;
}

Eigen::Matrix2d PlasticElasticMethod::calcDeviatorTensor(
    const std::size_t cell, const std::size_t quad, const std::size_t material,
    const double rhoLocal, const double jacobianDet,
    const Eigen::Matrix2d &jacobianInv, const Eigen::VectorXd &uCalc,
    const Eigen::VectorXd &deviatorStressTensor_00_calc,
    const Eigen::VectorXd &deviatorStressTensor_01_calc,
    const Eigen::VectorXd &deviatorStressTensor_10_calc,
    const Eigen::VectorXd &deviatorStressTensor_11_calc) {
  Eigen::Matrix2d deviator_stress_tensor;
  deviator_stress_tensor.setZero();

  for (std::size_t node = 0; node < kNumberOfThermodynamicPointsPerCell;
       node++) {
    const std::size_t thermo_index =
        getThermodynamicIndexFromCell(cell, node, material);

    const double s00_local = deviatorStressTensor_00_calc(thermo_index);
    const double s01_local = deviatorStressTensor_01_calc(thermo_index);
    const double s10_local = deviatorStressTensor_10_calc(thermo_index);
    const double s11_local = deviatorStressTensor_11_calc(thermo_index);

    const double node_basis = thermodynamicBasisQuadValues(node, quad);

    deviator_stress_tensor(0, 0) += s00_local * node_basis;
    deviator_stress_tensor(0, 1) += s01_local * node_basis;
    deviator_stress_tensor(1, 0) += s10_local * node_basis;
    deviator_stress_tensor(1, 1) += s11_local * node_basis;
  }

  Eigen::Matrix2d velocityGrad;
  velocityGrad.setZero();

  for (std::size_t kinematicNode = 0;
       kinematicNode < kNumberOfKinematicPointsPerCell; kinematicNode++) {
    const double kinematicBasisdx =
        kinematicBasisdxQuadValues(kinematicNode, quad);
    const double kinematicBasisdy =
        kinematicBasisdyQuadValues(kinematicNode, quad);

    for (std::size_t direction = 0; direction < kSolverDimention; direction++) {
      const std::size_t kinematicIndex =
          getKinematicIndexFromCell(cell, kinematicNode, direction);
      const double uNode = uCalc(kinematicIndex);
      velocityGrad(0, direction) += uNode * kinematicBasisdx;
      velocityGrad(1, direction) += uNode * kinematicBasisdy;
    }
  }
  velocityGrad = jacobianInv * velocityGrad;
  // velocityGrad.transposeInPlace();
  Eigen::Matrix2d symVelocityGrad = (velocityGrad + velocityGrad.transpose());
  Eigen::Matrix2d spinTensor = (velocityGrad - velocityGrad.transpose());
  const double velocityScalarGrad = velocityGrad(0, 0) + velocityGrad(1, 1);

  const double shear_modulus = problem.plastic_params->shear_modulus;

  const Eigen::Matrix2d term1 =
      shear_modulus * (symVelocityGrad - (2.0 / 3.0) * velocityScalarGrad *
                                             Eigen::Matrix2d::Identity());
  const Eigen::Matrix2d term2 = -0.5 * (deviator_stress_tensor * spinTensor -
                                        spinTensor * deviator_stress_tensor);

  const Eigen::Matrix2d g_local = term1 + term2;

  for (std::size_t thermo_node = 0;
       thermo_node < kNumberOfThermodynamicPointsPerCell; thermo_node++) {
    const double quad_weight = quadWeights[quad];

    const double g00_local = g_local(0, 0);
    const double g01_local = g_local(0, 1);
    const double g10_local = g_local(1, 0);
    const double g11_local = g_local(1, 1);

    const double thermo_basis_value =
        thermodynamicBasisQuadValues(thermo_node, quad);

    const double common =
        quad_weight * rhoLocal * thermo_basis_value * jacobianDet;

    const double g00_value = common * g00_local;
    const double g01_value = common * g01_local;
    const double g10_value = common * g10_local;
    const double g11_value = common * g11_local;

    const std::size_t thermo_index =
        getThermodynamicIndexFromCell(cell, thermo_node, material);

    deviator_rate_00(thermo_index) += g00_value;
    deviator_rate_01(thermo_index) += g01_value;
    deviator_rate_10(thermo_index) += g10_value;
    deviator_rate_11(thermo_index) += g11_value;
  }

  return deviator_stress_tensor;
}

void PlasticElasticMethod::resolveBoundaryKinematicMassMatrix() {
  resolveLeftBoundaryKinematicMassMatrix();
  resolveTopBoundaryKinematicMassMatrix();
  resolveRightBoundaryKinematicMassMatrix();
  resolveBottomBoundaryKinematicMassMatrix();
}
void PlasticElasticMethod::resolveLeftBoundaryKinematicMassMatrix() {
  switch (problem.leftBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (std::size_t cell = 0; cell < yCells; cell++) {
        for (std::size_t nodei = 0;
             nodei < kNumberOfKinematicPointsPerCellPerDimention; nodei++) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            constexpr std::size_t direction = 0;
            const std::size_t nodeiIndex =
                getKinematicIndexFromCell(cell, nodei, direction);
            const std::size_t nodejIndex =
                getKinematicIndexFromCell(cell, nodej, direction);
            kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) = 0.0;
            kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (std::size_t cell = 0; cell < yCells; cell++) {
        for (std::size_t nodei = 0;
             nodei < kNumberOfKinematicPointsPerCellPerDimention; nodei++) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            for (std::size_t direction = 0; direction < kSolverDimention;
                 direction++) {
              const std::size_t nodeiIndex =
                  getKinematicIndexFromCell(cell, nodei, direction);
              const std::size_t nodejIndex =
                  getKinematicIndexFromCell(cell, nodej, direction);
              kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) = 0.0;
              kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) = 0.0;
            }
          }
        }
      }
      break;
    case BoundaryType::eTransmissive:
      for (std::size_t cell = 0; cell < yCells; cell++) {
        for (std::size_t nodei = 0;
             nodei < kNumberOfKinematicPointsPerCellPerDimention; nodei++) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            for (std::size_t direction = 0; direction < kSolverDimention;
                 direction++) {
              const std::size_t nodeiIndex =
                  getKinematicIndexFromCell(cell, nodei, direction);
              const std::size_t nodejIndex =
                  getKinematicIndexFromCell(cell, nodej, direction);
              kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) *= 2.0;
              kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) *= 2.0;
            }
          }
        }
      }
      break;
  }
}
void PlasticElasticMethod::resolveTopBoundaryKinematicMassMatrix() {
  switch (problem.topBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (std::size_t cell = yCells - 1; cell < kNumberOfCells;
           cell += yCells) {
        for (std::size_t nodei =
                 kNumberOfKinematicPointsPerCellPerDimention - 1;
             nodei < kNumberOfKinematicPointsPerCell;
             nodei += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            constexpr std::size_t direction = 1;
            const std::size_t nodeiIndex =
                getKinematicIndexFromCell(cell, nodei, direction);
            const std::size_t nodejIndex =
                getKinematicIndexFromCell(cell, nodej, direction);
            kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) = 0.0;
            kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (std::size_t cell = yCells - 1; cell < kNumberOfCells;
           cell += yCells) {
        for (std::size_t nodei =
                 kNumberOfKinematicPointsPerCellPerDimention - 1;
             nodei < kNumberOfKinematicPointsPerCell;
             nodei += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            for (std::size_t direction = 0; direction < kSolverDimention;
                 direction++) {
              const std::size_t nodeiIndex =
                  getKinematicIndexFromCell(cell, nodei, direction);
              const std::size_t nodejIndex =
                  getKinematicIndexFromCell(cell, nodej, direction);
              kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) = 0.0;
              kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) = 0.0;
            }
          }
        }
      }
      break;
    case BoundaryType::eTransmissive:
      for (std::size_t cell = yCells - 1; cell < kNumberOfCells;
           cell += yCells) {
        for (std::size_t nodei =
                 kNumberOfKinematicPointsPerCellPerDimention - 1;
             nodei < kNumberOfKinematicPointsPerCell;
             nodei += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            for (std::size_t direction = 0; direction < kSolverDimention;
                 direction++) {
              const std::size_t nodeiIndex =
                  getKinematicIndexFromCell(cell, nodei, direction);
              const std::size_t nodejIndex =
                  getKinematicIndexFromCell(cell, nodej, direction);
              kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) *= 2.0;
              kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) *= 2.0;
            }
          }
        }
      }
      break;
  }
}
void PlasticElasticMethod::resolveRightBoundaryKinematicMassMatrix() {
  switch (problem.rightBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (std::size_t cell = (xCells - 1) * yCells; cell < kNumberOfCells;
           cell++) {
        for (std::size_t nodei =
                 (kNumberOfKinematicPointsPerCellPerDimention - 1) *
                 kNumberOfKinematicPointsPerCellPerDimention;
             nodei < kNumberOfKinematicPointsPerCell; nodei++) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            constexpr std::size_t direction = 0;
            const std::size_t nodeiIndex =
                getKinematicIndexFromCell(cell, nodei, direction);
            const std::size_t nodejIndex =
                getKinematicIndexFromCell(cell, nodej, direction);
            kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) = 0.0;
            kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (std::size_t cell = (xCells - 1) * yCells; cell < kNumberOfCells;
           cell++) {
        for (std::size_t nodei =
                 (kNumberOfKinematicPointsPerCellPerDimention - 1) *
                 kNumberOfKinematicPointsPerCellPerDimention;
             nodei < kNumberOfKinematicPointsPerCell; nodei++) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            for (std::size_t direction = 0; direction < kSolverDimention;
                 direction++) {
              const std::size_t nodeiIndex =
                  getKinematicIndexFromCell(cell, nodei, direction);
              const std::size_t nodejIndex =
                  getKinematicIndexFromCell(cell, nodej, direction);
              kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) = 0.0;
              kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) = 0.0;
            }
          }
        }
      }
      break;
    case BoundaryType::eTransmissive:
      for (std::size_t cell = (xCells - 1) * yCells; cell < kNumberOfCells;
           cell++) {
        for (std::size_t nodei =
                 (kNumberOfKinematicPointsPerCellPerDimention - 1) *
                 kNumberOfKinematicPointsPerCellPerDimention;
             nodei < kNumberOfKinematicPointsPerCell; nodei++) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            for (std::size_t direction = 0; direction < kSolverDimention;
                 direction++) {
              const std::size_t nodeiIndex =
                  getKinematicIndexFromCell(cell, nodei, direction);
              const std::size_t nodejIndex =
                  getKinematicIndexFromCell(cell, nodej, direction);
              kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) *= 2.0;
              kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) *= 2.0;
            }
          }
        }
      }
      break;
  }
}
void PlasticElasticMethod::resolveBottomBoundaryKinematicMassMatrix() {
  switch (problem.bottomBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (std::size_t cell = 0; cell < kNumberOfCells; cell += yCells) {
        for (std::size_t nodei = 0; nodei < kNumberOfKinematicPointsPerCell;
             nodei += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            constexpr std::size_t direction = 1;
            const std::size_t nodeiIndex =
                getKinematicIndexFromCell(cell, nodei, direction);
            const std::size_t nodejIndex =
                getKinematicIndexFromCell(cell, nodej, direction);

            kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) = 0.0;
            kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (std::size_t cell = 0; cell < kNumberOfCells; cell += yCells) {
        for (std::size_t nodei = 0; nodei < kNumberOfKinematicPointsPerCell;
             nodei += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            for (std::size_t direction = 0; direction < kSolverDimention;
                 direction++) {
              const std::size_t nodeiIndex =
                  getKinematicIndexFromCell(cell, nodei, direction);
              const std::size_t nodejIndex =
                  getKinematicIndexFromCell(cell, nodej, direction);

              kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) = 0.0;
              kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) = 0.0;
            }
          }
        }
      }
      break;
    case BoundaryType::eTransmissive:
      for (std::size_t cell = 0; cell < kNumberOfCells; cell += yCells) {
        for (std::size_t nodei = 0; nodei < kNumberOfKinematicPointsPerCell;
             nodei += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t nodej = 0; nodej < kNumberOfKinematicPointsPerCell;
               nodej++) {
            for (std::size_t direction = 0; direction < kSolverDimention;
                 direction++) {
              const std::size_t nodeiIndex =
                  getKinematicIndexFromCell(cell, nodei, direction);
              const std::size_t nodejIndex =
                  getKinematicIndexFromCell(cell, nodej, direction);

              kinematicMassMatrix.coeffRef(nodeiIndex, nodejIndex) *= 2.0;
              kinematicMassMatrix.coeffRef(nodejIndex, nodeiIndex) *= 2.0;
            }
          }
        }
      }
      break;
  }
}

void PlasticElasticMethod::resolveBoundaryForceVector(Eigen::VectorXd &FuCalc) {
  resolveLeftBoundaryForceVector(FuCalc);
  resolveTopBoundaryForceVector(FuCalc);
  resolveRightBoundaryForceVector(FuCalc);
  resolveBottomBoundaryForceVector(FuCalc);
}
void PlasticElasticMethod::resolveLeftBoundaryForceVector(
    Eigen::VectorXd &FuCalc) {
  switch (problem.leftBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (std::size_t cell = 0; cell < yCells; cell++) {
        for (std::size_t node = 0;
             node < kNumberOfKinematicPointsPerCellPerDimention; node++) {
          constexpr std::size_t direction = 0;
          const std::size_t kinematicIndex =
              getKinematicIndexFromCell(cell, node, direction);
          FuCalc(kinematicIndex) = 0.0;
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (std::size_t cell = 0; cell < yCells; cell++) {
        for (std::size_t node = 0;
             node < kNumberOfKinematicPointsPerCellPerDimention; node++) {
          for (std::size_t direction = 0; direction < kSolverDimention;
               direction++) {
            const std::size_t kinematicIndex =
                getKinematicIndexFromCell(cell, node, direction);
            FuCalc(kinematicIndex) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eTransmissive:
      for (std::size_t cell = 0; cell < yCells; cell++) {
        for (std::size_t node = 0;
             node < kNumberOfKinematicPointsPerCellPerDimention; node++) {
          for (std::size_t direction = 0; direction < kSolverDimention;
               direction++) {
            const std::size_t kinematicIndex =
                getKinematicIndexFromCell(cell, node, direction);
            FuCalc(kinematicIndex) = 0.0;
          }
        }
      }
      break;
  }
}
void PlasticElasticMethod::resolveTopBoundaryForceVector(
    Eigen::VectorXd &FuCalc) {
  assert(yCells >= 1);
  switch (problem.topBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (std::size_t cell = yCells - 1; cell < kNumberOfCells;
           cell += yCells) {
        for (std::size_t node = kNumberOfKinematicPointsPerCellPerDimention - 1;
             node < kNumberOfKinematicPointsPerCell;
             node += kNumberOfKinematicPointsPerCellPerDimention) {
          constexpr std::size_t direction = 1;
          const std::size_t kinematicIndex =
              getKinematicIndexFromCell(cell, node, direction);
          FuCalc(kinematicIndex) = 0.0;
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (std::size_t cell = yCells - 1; cell < kNumberOfCells;
           cell += yCells) {
        for (std::size_t node = kNumberOfKinematicPointsPerCellPerDimention - 1;
             node < kNumberOfKinematicPointsPerCell;
             node += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t direction = 0; direction < kSolverDimention;
               direction++) {
            const std::size_t kinematicIndex =
                getKinematicIndexFromCell(cell, node, direction);
            FuCalc(kinematicIndex) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eTransmissive:
      for (std::size_t cell = yCells - 1; cell < kNumberOfCells;
           cell += yCells) {
        for (std::size_t node = kNumberOfKinematicPointsPerCellPerDimention - 1;
             node < kNumberOfKinematicPointsPerCell;
             node += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t direction = 0; direction < kSolverDimention;
               direction++) {
            const std::size_t kinematicIndex =
                getKinematicIndexFromCell(cell, node, direction);
            FuCalc(kinematicIndex) = 0.0;
          }
        }
      }
      break;
  }
}
void PlasticElasticMethod::resolveRightBoundaryForceVector(
    Eigen::VectorXd &FuCalc) {
  assert(xCells >= 1);
  assert(kNumberOfKinematicPointsPerCellPerDimention >= 1);
  switch (problem.rightBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (std::size_t cell = (xCells - 1) * yCells; cell < kNumberOfCells;
           cell++) {
        for (std::size_t node =
                 (kNumberOfKinematicPointsPerCellPerDimention - 1) *
                 kNumberOfKinematicPointsPerCellPerDimention;
             node < kNumberOfKinematicPointsPerCell; node++) {
          constexpr std::size_t direction = 0;
          const std::size_t kinematicIndex =
              getKinematicIndexFromCell(cell, node, direction);
          FuCalc(kinematicIndex) = 0.0;
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (std::size_t cell = (xCells - 1) * yCells; cell < kNumberOfCells;
           cell++) {
        for (std::size_t node =
                 (kNumberOfKinematicPointsPerCellPerDimention - 1) *
                 kNumberOfKinematicPointsPerCellPerDimention;
             node < kNumberOfKinematicPointsPerCell; node++) {
          for (std::size_t direction = 0; direction < kSolverDimention;
               direction++) {
            const std::size_t kinematicIndex =
                getKinematicIndexFromCell(cell, node, direction);
            FuCalc(kinematicIndex) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eTransmissive:
      for (std::size_t cell = (xCells - 1) * yCells; cell < kNumberOfCells;
           cell++) {
        for (std::size_t node =
                 (kNumberOfKinematicPointsPerCellPerDimention - 1) *
                 kNumberOfKinematicPointsPerCellPerDimention;
             node < kNumberOfKinematicPointsPerCell; node++) {
          for (std::size_t direction = 0; direction < kSolverDimention;
               direction++) {
            const std::size_t kinematicIndex =
                getKinematicIndexFromCell(cell, node, direction);
            FuCalc(kinematicIndex) = 0.0;
          }
        }
      }
      break;
  }
}
void PlasticElasticMethod::resolveBottomBoundaryForceVector(
    Eigen::VectorXd &FuCalc) {
  switch (problem.bottomBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (std::size_t cell = 0; cell < kNumberOfCells; cell += yCells) {
        for (std::size_t node = 0; node < kNumberOfKinematicPointsPerCell;
             node += kNumberOfKinematicPointsPerCellPerDimention) {
          constexpr std::size_t direction = 1;
          const std::size_t kinematicIndex =
              getKinematicIndexFromCell(cell, node, direction);
          FuCalc(kinematicIndex) = 0.0;
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (std::size_t cell = 0; cell < kNumberOfCells; cell += yCells) {
        for (std::size_t node = 0; node < kNumberOfKinematicPointsPerCell;
             node += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t direction = 0; direction < kSolverDimention;
               direction++) {
            const std::size_t kinematicIndex =
                getKinematicIndexFromCell(cell, node, direction);
            FuCalc(kinematicIndex) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eTransmissive:
      for (std::size_t cell = 0; cell < kNumberOfCells; cell += yCells) {
        for (std::size_t node = 0; node < kNumberOfKinematicPointsPerCell;
             node += kNumberOfKinematicPointsPerCellPerDimention) {
          for (std::size_t direction = 0; direction < kSolverDimention;
               direction++) {
            const std::size_t kinematicIndex =
                getKinematicIndexFromCell(cell, node, direction);
            FuCalc(kinematicIndex) = 0.0;
          }
        }
      }
      break;
  }
}

void PlasticElasticMethod::RK2Step() {
  while (true) {
    tau = std::numeric_limits<double>::max();
    calcForceMatrix(x, u, volFrac, e, stress_deviator_00, stress_deviator_01,
                    stress_deviator_10, stress_deviator_11);
    if (dt >= tau) {
      dt = beta1 * tau;
    }
    stress_deviator_00_05 = thermoMassMatrixInv * deviator_rate_00;
    stress_deviator_00_05 *= 0.5 * dt;
    stress_deviator_00_05 += stress_deviator_00;

    stress_deviator_01_05 = thermoMassMatrixInv * deviator_rate_01;
    stress_deviator_01_05 *= 0.5 * dt;
    stress_deviator_01_05 += stress_deviator_01;

    stress_deviator_10_05 = thermoMassMatrixInv * deviator_rate_10;
    stress_deviator_10_05 *= 0.5 * dt;
    stress_deviator_10_05 += stress_deviator_10;

    stress_deviator_11_05 = thermoMassMatrixInv * deviator_rate_11;
    stress_deviator_11_05 *= 0.5 * dt;
    stress_deviator_11_05 += stress_deviator_11;

    equivalent_plastic_strain_05 = equivalent_plastic_strain;
    radialReturn(stress_deviator_00_05, stress_deviator_01_05,
                 stress_deviator_10_05, stress_deviator_11_05,
                 equivalent_plastic_strain_05);

    Fu = forceMatrix * Eigen::VectorXd::Ones(kNumberOfMaterials *
                                             kNumberOfThermodynamicPointsTotal);
    resolveBoundaryForceVector(Fu);

    u05 = kinematicMassMatrixSolver.solve(Fu);
    if (kinematicMassMatrixSolver.info() != Eigen::Success) {
      throw std::runtime_error("failed to solve kinematic mass matrix");
    }

    u05 *= -0.5 * dt;
    u05 += u;

    Fe = forceMatrix.transpose() * u05 - energyExchangeRate;
    e05 = thermoMassMatrixInv * Fe;
    e05 *= 0.5 * dt;
    e05 += e;

    x05 = x + 0.5 * dt * u05;

    volFrac05 = volFrac + 0.5 * dt * volFracRate;

    calcForceMatrix(x05, u05, volFrac05, e05, stress_deviator_00_05,
                    stress_deviator_01_05, stress_deviator_10_05,
                    stress_deviator_11_05);
    if (dt >= tau) {
      dt = beta1 * tau;
      continue;
    } else {
      break;
    }
  }
  stress_deviator_00_05 = thermoMassMatrixInv * deviator_rate_00;
  stress_deviator_00_05 *= dt;
  stress_deviator_00_05 += stress_deviator_00;
  stress_deviator_00.swap(stress_deviator_00_05);

  stress_deviator_01_05 = thermoMassMatrixInv * deviator_rate_01;
  stress_deviator_01_05 *= dt;
  stress_deviator_01_05 += stress_deviator_01;
  stress_deviator_01.swap(stress_deviator_01_05);

  stress_deviator_10_05 = thermoMassMatrixInv * deviator_rate_10;
  stress_deviator_10_05 *= dt;
  stress_deviator_10_05 += stress_deviator_10;
  stress_deviator_10.swap(stress_deviator_10_05);

  stress_deviator_11_05 = thermoMassMatrixInv * deviator_rate_11;
  stress_deviator_11_05 *= dt;
  stress_deviator_11_05 += stress_deviator_11;
  stress_deviator_11.swap(stress_deviator_11_05);

  equivalent_plastic_strain.swap(equivalent_plastic_strain_05);
  radialReturn(stress_deviator_00, stress_deviator_01, stress_deviator_10,
               stress_deviator_11, equivalent_plastic_strain);

  Fu = forceMatrix * Eigen::VectorXd::Ones(kNumberOfMaterials *
                                           kNumberOfThermodynamicPointsTotal);
  resolveBoundaryForceVector(Fu);

  u05 = kinematicMassMatrixSolver.solve(Fu);
  if (kinematicMassMatrixSolver.info() != Eigen::Success) {
    throw std::runtime_error("failed to solve kinematic mass matrix");
  }

  u05 *= -dt;
  u05 += u;
  u.swap(u05);
  u05 += u;
  u05 *= 0.5;

  Fe = forceMatrix.transpose() * u05 - energyExchangeRate;
  e05 = thermoMassMatrixInv * Fe;
  e05 *= dt;
  e05 += e;

  x05 = x + dt * u05;

  volFrac05 = volFrac + dt * volFracRate;

  x.swap(x05);
  volFrac.swap(volFrac05);
  e.swap(e05);

  t += dt;
  if (dt <= gamma * tau) {
    dt = beta2 * dt;
  }
}

void PlasticElasticMethod::dumpPlasticityThreshold(
    const std::filesystem::path &dataOutputDir) const {
  for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
    const std::string filenamePlasticityThreshold =
        "plasticity_threshold" + std::to_string(material) + ".fem";

    std::ofstream ofsPlasticityThreshold(
        dataOutputDir / filenamePlasticityThreshold, std::ios::binary);

    for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
      for (std::size_t out = 0; out < kNumberOfOutputPointsPerCell; out++) {
        double s00_out = 0.0;
        double s01_out = 0.0;
        double s10_out = 0.0;
        double s11_out = 0.0;
        double s22_out = 0.0;

        for (std::size_t node = 0; node < kNumberOfThermodynamicPointsPerCell;
             node++) {
          const std::size_t nodeIndex =
              getThermodynamicIndexFromCell(cell, node, material);

          const double s00_local = stress_deviator_00(nodeIndex);
          const double s01_local = stress_deviator_01(nodeIndex);
          const double s10_local = stress_deviator_10(nodeIndex);
          const double s11_local = stress_deviator_11(nodeIndex);
          const double s22_local = -s00_local - s11_local;

          const double nodeBasis = thermodynamicBasisOutputValues(node, out);

          s00_out += s00_local * nodeBasis;
          s01_out += s01_local * nodeBasis;
          s10_out += s10_local * nodeBasis;
          s11_out += s11_local * nodeBasis;
          s22_out += s22_local * nodeBasis;
        }
        const double stress_measure_out = std::sqrt(
            1.5 * (s00_out * s00_out + s01_out * s01_out + s10_out * s10_out +
                   s11_out * s11_out + s22_out * s22_out));

        const double equivalent_stress =
            problem.plastic_params->equivalent_stress;
        const double plasticity_threshold =
            stress_measure_out / equivalent_stress;

        ofsPlasticityThreshold.write(
            reinterpret_cast<const char *>(&plasticity_threshold),
            sizeof(double));
      }
    }
  }
}

void PlasticElasticMethod::dumpDeviatorStressTensor(
    const std::filesystem::path &dataOutputDir) const {
  for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
    const std::string filenameS00 = "s00_" + std::to_string(material) + ".fem";
    const std::string filenameS01 = "s01_" + std::to_string(material) + ".fem";
    const std::string filenameS10 = "s10_" + std::to_string(material) + ".fem";
    const std::string filenameS11 = "s11_" + std::to_string(material) + ".fem";

    std::ofstream ofsS00(dataOutputDir / filenameS00, std::ios::binary);
    std::ofstream ofsS01(dataOutputDir / filenameS01, std::ios::binary);
    std::ofstream ofsS10(dataOutputDir / filenameS10, std::ios::binary);
    std::ofstream ofsS11(dataOutputDir / filenameS11, std::ios::binary);

    for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
      for (std::size_t out = 0; out < kNumberOfOutputPointsPerCell; out++) {
        double s00_out = 0.0;
        double s01_out = 0.0;
        double s10_out = 0.0;
        double s11_out = 0.0;

        for (std::size_t node = 0; node < kNumberOfThermodynamicPointsPerCell;
             node++) {
          const std::size_t nodeIndex =
              getThermodynamicIndexFromCell(cell, node, material);

          const double s00_local = stress_deviator_00(nodeIndex);
          const double s01_local = stress_deviator_01(nodeIndex);
          const double s10_local = stress_deviator_10(nodeIndex);
          const double s11_local = stress_deviator_11(nodeIndex);

          const double nodeBasis = thermodynamicBasisOutputValues(node, out);

          s00_out += s00_local * nodeBasis;
          s01_out += s01_local * nodeBasis;
          s10_out += s10_local * nodeBasis;
          s11_out += s11_local * nodeBasis;
        }

        ofsS00.write(reinterpret_cast<const char *>(&s00_out), sizeof(double));
        ofsS01.write(reinterpret_cast<const char *>(&s01_out), sizeof(double));
        ofsS10.write(reinterpret_cast<const char *>(&s10_out), sizeof(double));
        ofsS11.write(reinterpret_cast<const char *>(&s11_out), sizeof(double));
      }
    }
  }
}

void PlasticElasticMethod::dumpEquivalentPlasticStrain(
    const std::filesystem::path &dataOutputDir) const {
  for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
    const std::string filenameEPS = "EPS" + std::to_string(material) + ".fem";

    std::ofstream ofsEPS(dataOutputDir / filenameEPS, std::ios::binary);

    for (std::size_t cell = 0; cell < kNumberOfCells; cell++) {
      for (std::size_t out = 0; out < kNumberOfOutputPointsPerCell; out++) {
        double EPS_out = 0.0;

        for (std::size_t thermo_node = 0;
             thermo_node < kNumberOfThermodynamicPointsPerCell; thermo_node++) {
          const std::size_t thermo_index =
              getThermodynamicIndexFromCell(cell, thermo_node, material);

          const double EPS_local = equivalent_plastic_strain(thermo_index);

          const double thermo_basis =
              thermodynamicBasisOutputValues(thermo_node, out);

          EPS_out += EPS_local * thermo_basis;
        }

        ofsEPS.write(reinterpret_cast<const char *>(&EPS_out), sizeof(double));
      }
    }
  }
}

Eigen::Matrix2d PlasticElasticMethod::calcJacobian(
    const std::size_t cell, const std::size_t quad,
    const Eigen::VectorXd &xCalc) const {
  Eigen::Matrix2d output;
  output.setZero();

  for (std::size_t node = 0; node < kNumberOfKinematicPointsPerCell; node++) {
    for (std::size_t direction = 0; direction < kSolverDimention; direction++) {
      const std::size_t nodeIndex =
          getKinematicIndexFromCell(cell, node, direction);
      const double xLocal = xCalc(nodeIndex);
      const double kinematicBasisdx = kinematicBasisdxQuadValues(node, quad);
      const double kinematicBasisdy = kinematicBasisdyQuadValues(node, quad);
      output(0, direction) += xLocal * kinematicBasisdx;
      output(1, direction) += xLocal * kinematicBasisdy;
    }
  }

  return output;
}
