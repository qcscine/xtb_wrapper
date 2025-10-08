/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <stdlib.h>
#include <algorithm>
#include <boost/dll/runtime_symbol_info.hpp>
#include <vector>

using namespace testing;
namespace Scine {
namespace Xtb {
namespace Tests {

class CalculatorTest : public TestWithParam<std::string> {
 public:
  static const std::vector<std::string> methods;
  static const std::map<std::string, double> energyMap;
  static const std::map<std::string, Eigen::MatrixXd> gradientMap;
  std::shared_ptr<Core::Calculator> calculator;
  boost::filesystem::path pathToResource;
  Utils::AtomCollection molecule;
  std::string method;

  static void SetUpTestSuite() {
    auto xtbPath = boost::dll::program_location().parent_path() / "../../xtb-src";
    // This is necessary for the param_gfn0-xtb.txt to be found, ensuring the GFN0 test passes
    setenv("XTBPATH", xtbPath.string().c_str(), 1);
  }

  static void TearDownTestSuite() {
    unsetenv("XTBPATH");
  }

 private:
  void SetUp() final {
    pathToResource = boost::dll::program_location().parent_path() / "Resources";
    auto xtbPath = boost::dll::program_location().parent_path() / "../../xtb-src";
    loadAtomCollection(molecule, "h2o.xyz");
    method = GetParam();
    auto& manager = Core::ModuleManager::getInstance();
    calculator = manager.get<Core::Calculator>(method);
    calculator->settings().modifyInt(Utils::SettingsNames::externalProgramNProcs, 1);
  }

  void loadAtomCollection(Utils::AtomCollection& collection, const std::string& filename) {
    auto filePath = pathToResource / filename;
    std::ifstream file(filePath.string());
    if (file) {
      std::stringstream buffer;
      buffer << file.rdbuf();
      collection = Utils::XyzStreamHandler::read(buffer);
    }
    else {
      throw std::runtime_error("Failed to open file: " + filePath.string());
    }
  }
};

const std::vector<std::string> CalculatorTest::methods = {"gfn0", "gfn1", "gfn2", "gfnff"};
// generated with ./xtb --gfn [0,1,2] --grad h2o.xyz and ./xtb --gfnff --grad h2o.xyz
const std::map<std::string, double> CalculatorTest::energyMap = {
    {"gfn0", -4.366510048328}, {"gfn1", -5.768497364658}, {"gfn2", -5.070194077901}, {"gfnff", -0.326794799815}};
const std::map<std::string, Eigen::MatrixXd> CalculatorTest::gradientMap = {
    {"gfn0", (Eigen::MatrixXd(3, 3) << 5.4633883002136E-03, 7.4807641149800E-03, 6.0951910779116E-17, -1.2642072114675E-02,
              4.0419268248342E-03, -6.5691326914328E-17, 7.1786838144614E-03, -1.1522690939814E-02, 9.9602740682517E-18)
                 .finished()},

    {"gfn1", (Eigen::MatrixXd(3, 3) << 3.1619965443865E-03, 9.6161804877529E-03, -4.2718253409498E-17, -8.8498732169638E-03,
              -1.5799556317044E-03, 5.0020782932492E-17, 5.6878766725773E-03, -8.0362248560486E-03, -7.3025295229946E-18)
                 .finished()},

    {"gfn2", (Eigen::MatrixXd(3, 3) << 3.3577653820128E-03, -6.3999842264807E-03, 7.3351985598374E-17, -1.0037798452429E-02,
              1.5383790457195E-02, -5.9913623759769E-17, 6.6800330704163E-03, -8.9838062307138E-03, -1.3438361838606E-17)
                 .finished()},

    {"gfnff", (Eigen::MatrixXd(3, 3) << 1.2469287336730E-02, 1.6607504316817E-02, 0.0000000000000E+00, -2.0459059320711E-02,
               2.8322750107321E-03, 0.0000000000000E+00, 7.9897719839803E-03, -1.9439779327549E-02, 0.0000000000000E+00)
                  .finished()},
};

INSTANTIATE_TEST_SUITE_P(XtbWrapper, CalculatorTest, ValuesIn(CalculatorTest::methods));

/*
Please note that depending on the GCC version in use, this test may fail due to known upstream issues.
Refer to the following links for more information:
- https://github.com/grimme-lab/xtb/pull/907
- https://github.com/grimme-lab/xtb/pull/1292
- https://github.com/grimme-lab/xtb/issues/1326
*/
TEST_P(CalculatorTest, TestCalculations) {
  calculator->setStructure(molecule);
  calculator->settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 1);
  calculator->setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  auto results = calculator->calculate("");
  EXPECT_NEAR(results.get<Utils::Property::Energy>(), CalculatorTest::energyMap.at(method), 1e-6)
      << "Energy does not match expected value!";
  ASSERT_TRUE(results.get<Utils::Property::Gradients>().isApprox(CalculatorTest::gradientMap.at(method), 1e-5))
      << "Gradient does not match expected value!";
}

TEST_P(CalculatorTest, CheckResultsClearing1) {
  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");
  auto structure = Utils::XyzStreamHandler::read(stream);
  calculator->results().set<Utils::Property::Energy>(42.0);
  calculator->setStructure(structure);
  ASSERT_FALSE(calculator->results().has<Utils::Property::Energy>());
}

TEST_P(CalculatorTest, CloneInterfaceWorksCorrectly) {
  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");
  auto structure = Utils::XyzStreamHandler::read(stream);
  calculator->settings().modifyInt(Utils::SettingsNames::externalProgramNProcs, 2);

  calculator->setStructure(structure);
  calculator->results().set<Utils::Property::Energy>(42.0);
  auto newCalculator = std::shared_ptr<Core::Calculator>(calculator->clone());

  ASSERT_THAT(calculator->getPositions()(3, 1), Eq(newCalculator->getPositions()(3, 1)));
  ASSERT_THAT(calculator->getPositions()(4, 2), Eq(newCalculator->getPositions()(4, 2)));
  ASSERT_THAT(calculator->results().get<Utils::Property::Energy>(),
              Eq(newCalculator->results().get<Utils::Property::Energy>()));
  ASSERT_THAT(calculator->settings().getInt(Utils::SettingsNames::externalProgramNProcs),
              Eq(newCalculator->settings().getInt(Utils::SettingsNames::externalProgramNProcs)));
}

} // namespace Tests
} // namespace Xtb
} // namespace Scine
