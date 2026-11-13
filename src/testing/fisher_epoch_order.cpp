#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gsl/gsl_deriv.h>

#include "../fisher.cpp"

namespace {

void require(bool condition, const std::string& message)
{
  if (!condition)
    {
      std::cerr << "FAIL: " << message << '\n';
      std::exit(1);
    }
}

bool close(double actual, double expected)
{
  // Allow tiny solver roundoff while remaining well below the stale-state shift.
  return std::isfinite(actual) && std::isfinite(expected) &&
    std::abs(actual - expected) <= 1e-7 +
    1e-6 * std::max(std::abs(actual), std::abs(expected));
}

struct QuietOutput
{
  std::ostringstream sink;
  std::streambuf* previous;
  QuietOutput() : previous(std::cout.rdbuf(sink.rdbuf())) {}
  ~QuietOutput() { std::cout.rdbuf(previous); }
};

struct OrderFixture
{
  static constexpr int observers = 2;
  static constexpr int epochsPerObserver = 4;

  filekeywords parameters{};
  event data{};
  obsfilekeywords world[observers]{};
  VBMicrolensing solver;

  OrderFixture()
  {
    parameters.pllxMultiplyer = 1.0;
    parameters.LD_GAMMA = 0.0;

    data.id = 1;
    data.t0 = 92.0;
    data.tref = 90.0;
    data.tE_r = 28.0;
    data.u0 = 0.31;
    data.thE = data.rE = 1.0;
    data.rs = 0.001;
    data.alpha = 27.0;
    data.piEN = 0.13;
    data.piEE = -0.07;
    data.piE = std::hypot(data.piEN, data.piEE);
    data.nepochs = observers * epochsPerObserver;
    data.vbm = &solver;
    solver.Tol = solver.RelTol = 1e-10;

    constexpr double q = 0.05;
    constexpr double separation = 1.4;
    constexpr double rotation = 0.43;
    data.p_q = {q};
    data.p_s0 = {separation};
    data.p_x0 = {separation * std::cos(rotation)};
    data.p_y0 = {separation * std::sin(rotation)};
    data.pllx.resize(observers);

    for (int obs = 0; obs < observers; ++obs)
      {
        // Each observer has its own local shift array and distinct track values.
        data.nepochsvec[obs] = obs * epochsPerObserver;
        parallax& track = data.pllx[obs];
        track.set_lb(0.02, -0.03);
        track.status = parallax::allready;
        track.vref = {0.0, 0.0, 0.0};
        track.aref = {0.001, 0.002, 0.0};
        track.provide_observables_NE(data.piEN, data.piEE,
                                     data.tE_r, data.thE);

        for (int point = 0; point < epochsPerObserver; ++point)
          {
            const double time = 58.0 + 19.0 * point + 0.37 * obs;
            const double dt = (time - data.tref) / 100.0;
            data.epoch.push_back(time);
            data.obsidx.push_back(obs);
            track.epochs.push_back(time + 2458849.0);
            track.Nshift.push_back(0.24 * dt * dt + 0.11 * obs);
            track.Eshift.push_back(0.13 * dt * dt * dt - 0.08 * obs);
            track.tshift.push_back(0.0);
            track.ushift.push_back(0.0);
          }
        track.compute_tushifts();
      }
    data.nepochsvec[observers] = data.nepochs;
  }
};

double rawT0Derivative(OrderFixture& fixture, int index)
{
  // Use Fisher's raw derivative path only as a repeatability probe; the test
  // makes no claim about derivative accuracy or a preferred numerical step.
  DerivContext context{0, index, &fixture.parameters, &fixture.data,
                       fixture.world};
  gsl_function function;
  function.function = deriv_wrapper;
  function.params = &context;
  double result = 0.0, error = 0.0;
  const int status = gsl_deriv_central(&function, 0.0, 1e-4,
                                       &result, &error);
  require(status == 0 && std::isfinite(result) && std::isfinite(error),
          "raw t0 derivative is finite for every fixture epoch");
  return result;
}

std::vector<double> baselineDerivatives()
{
  OrderFixture fixture;
  std::vector<double> result(fixture.data.nepochs);
  for (int index = 0; index < fixture.data.nepochs; ++index)
    result[index] = rawT0Derivative(fixture, index);
  return result;
}

std::vector<double> traverseWithParallaxPerturbations(bool reverse)
{
  OrderFixture fixture;
  std::vector<double> result(fixture.data.nepochs);
  for (int visit = 0; visit < fixture.data.nepochs; ++visit)
    {
      const int index = reverse ? fixture.data.nepochs - 1 - visit : visit;

      // Prime the same observer and epoch with a temporary piE change. The
      // next t0 derivative must rebuild its shift instead of reusing this one.
      // Alternate the two piE components so both direction and magnitude change.
      const int parallaxParameter = (index % 2 == 0) ? 7 : 8;
      DerivContext perturbation{parallaxParameter, index, &fixture.parameters,
                                &fixture.data, fixture.world};
      const double delta = (parallaxParameter == 7) ? 0.04 : -0.035;
      deriv_wrapper(delta, &perturbation);
      result[index] = rawT0Derivative(fixture, index);
    }
  return result;
}

void checkEpochTraversalOrder()
{
  const std::vector<double> baseline = baselineDerivatives();
  const std::vector<double> forward = traverseWithParallaxPerturbations(false);
  const std::vector<double> reverse = traverseWithParallaxPerturbations(true);

  for (std::size_t index = 0; index < baseline.size(); ++index)
    {
      if (!close(forward[index], baseline[index]) ||
          !close(reverse[index], baseline[index]))
        std::cerr << std::setprecision(17) << "epoch=" << index
                  << " baseline=" << baseline[index]
                  << " forward=" << forward[index]
                  << " reverse=" << reverse[index]
                  << " delta=" << (forward[index] - baseline[index]) << '\n';
      require(close(forward[index], baseline[index]),
              "forward traversal does not reuse a prior parallax shift");
      require(close(reverse[index], baseline[index]),
              "reverse traversal does not reuse a prior parallax shift");
      require(close(forward[index], reverse[index]),
              "the same raw epoch derivative is independent of traversal order");
    }
}

} // namespace

int main()
{
  {
    QuietOutput quiet;
    checkEpochTraversalOrder();
  }
  std::cout << "Fisher epoch-order check passed\n";
  return 0;
}
