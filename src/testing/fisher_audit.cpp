#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gsl/gsl_deriv.h>

#include "VBMicrolensingLibrary.h"
#include "lightcurveGenerator.h"

// This audit checks binary coordinate transforms, the general generator's
// CoM translation, and the fixed geometry used for Fisher flux terms.
static bool suppressFirstDerivative = false;
static int auditedCentral(const gsl_function*, double, double, double*, double*);

#define gsl_deriv_central auditedCentral
#include "../fisher.cpp"
#undef gsl_deriv_central

static int auditedCentral(const gsl_function* function, double x, double step,
                          double* result, double* error)
{
  const int status = gsl_deriv_central(function, x, step, result, error);
  const auto* context = static_cast<const DerivContext*>(function->params);
  // Keep the valid derivative; raise one reported error estimate to exercise the existing suppression path.
  if (suppressFirstDerivative && context->param == 0 && context->idx == 0)
    *error = 2.0 * std::abs(*result) + 1.0;
  return status;
}

namespace {

void require(bool condition, const std::string& message)
{
  if (!condition)
    {
      std::cerr << "FAIL: " << message << '\n';
      std::exit(1);
    }
}

bool close(double actual, double expected, double absoluteTolerance,
           double relativeTolerance = 0.0)
{
  return std::isfinite(actual) && std::isfinite(expected) &&
    std::abs(actual - expected) <= absoluteTolerance +
    relativeTolerance * std::max(std::abs(actual), std::abs(expected));
}

struct QuietOutput
{
  std::ostringstream sink;
  std::streambuf* previous;
  QuietOutput() : previous(std::cout.rdbuf(sink.rdbuf())) {}
  ~QuietOutput() { std::cout.rdbuf(previous); }
};

struct Fixture
{
  filekeywords param{};
  event data{};
  obsfilekeywords world[3];
  slcat sources, lenses;
  VBMicrolensing solver;

  explicit Fixture(double q = 0.05, bool withParallax = true,
                   bool liveLensOrbit = true)
  {
    param.numobservatories = 3;
    param.pllxMultiplyer = withParallax;
    param.LD_GAMMA = 0.0;
    param.simulation_zerotime = 2458849.0;
    param.num_lens_max = 2;
    param.lc_timeout = 0.0;
    param.vbm_tol = 1e-10;
    param.vbm_reltol = 1e-10;
    param.verbosity = 0;
    param.multiple_sources = 0;
    param.multiple_lenses = 0;
    param.lensOrbitalMotion = liveLensOrbit ? 1 : 0;

    data.id = 1;
    data.source = data.lens = 0;
    data.nsrc = 1;
    data.nlens = 2;
    data.nplanets = data.ncatalog_planets = 1;
    data.allsat = 0;
    data.numobservatories = 3;
    data.t0 = 100.0;
    data.tref = 95.0;
    data.tE_r = 25.0;
    data.u0 = 0.3;
    data.rE = data.thE = 1.0;
    data.rs = 0.001;
    data.ra = 4.6;
    data.dec = -0.5;
    data.nepochs = 33;
    data.piEN = 0.13;
    data.piEE = -0.07;
    data.piE = std::hypot(data.piEN, data.piEE);
    data.vbm = &solver;
    solver.Tol = solver.RelTol = 1e-10;

    data.p_mass = {q * 0.5};
    data.p_a = {1.4};
    data.p_e = {0.1};
    data.p_I = {43.0};
    data.p_L0 = {23.0};
    data.p_w = {17.0};
    data.p_O = {31.0};
    data.p_orbtype = {1.0};
    data.p_q = {q};
    for (auto* values : {&data.p_dL, &data.p_period, &data.p_s0,
                         &data.p_x0, &data.p_y0, &data.p_z0,
                         &data.p_dsdt, &data.p_dalphadt, &data.p_dxdt,
                         &data.p_dydt, &data.p_dzdt})
      values->assign(1, 0.0);

    data.pllx.resize(3);
    data.jdtimes.resize(3);
    data.Atrue.assign(data.nepochs, 0.0);
    data.Aerr.assign(data.nepochs, 0.003);
    data.mu_src.assign(1, std::vector<double>(data.nepochs));
    data.xsrc = data.ysrc = data.mu_src;
    data.xlens.assign(2, std::vector<double>(data.nepochs));
    data.ylens = data.xlens;
    data.tpeak.assign(1, 0.0);
    data.upeak.assign(1, 0.0);
    data.umin = 1e100;
    data.nosat.assign(data.nepochs, true);
    sources.DIST = 0;
    lenses.DIST = 0;
    lenses.MASS = 1;
    sources.data = {{8.0}};
    lenses.data = {{4.0, 0.5}};
    sources.mags = {std::vector<double>(32, 20.0)};
    lenses.mags = {std::vector<double>(32, 20.0)};

    constexpr int pointsPerObserver = 11;
    for (int obs = 0; obs < 3; ++obs)
      {
        world[obs].filter = 0;
        data.fs[obs] = 0.2 + 0.3 * obs;
        data.nepochsvec[obs] = obs * pointsPerObserver;
        auto& parallaxData = data.pllx[obs];
        parallaxData.set_lb(0.02, -0.03);
        parallaxData.status = parallax::allready;
        parallaxData.vref.assign(3, 0.0);
        parallaxData.aref = {0.001, 0.002, 0.0};
        parallaxData.provide_observables_NE(data.piEN, data.piEE,
                                             data.tE_r, data.thE);
        if (obs == 0)
          {
            data.murel_l = parallaxData.mul_h * parallaxData.murel_h;
            data.murel_b = parallaxData.mub_h * parallaxData.murel_h;
          }

        for (int point = 0; point < pointsPerObserver; ++point)
          {
            const double time = 20.0 + 16.0 * point + 0.31 * obs;
            const double dt = (time - data.tref) / 100.0;
            data.epoch.push_back(time);
            data.obsidx.push_back(obs);
            data.jdtimes[obs].push_back(time + param.simulation_zerotime);
            parallaxData.epochs.push_back(time + param.simulation_zerotime);
            parallaxData.Nshift.push_back(0.2 * dt * dt + 0.1 * obs);
            parallaxData.Eshift.push_back(0.1 * dt * dt * dt - 0.07 * obs);
            parallaxData.tshift.push_back(0.0);
            parallaxData.ushift.push_back(0.0);
          }
        parallaxData.compute_tushifts();
      }
    data.nepochsvec[3] = data.nepochs;

    std::ofstream noLog;
    QuietOutput quiet;
    lightcurveGenerator(&param, &data, world, &sources, &lenses, noLog);
    require(data.lcerror == 0, "general generator fixture completes");
  }
};

double explicitLensMagnification(const Fixture& fixture, int index)
{
  const event& ev = fixture.data;
  const int obs = ev.obsidx[index];
  const int local = index - ev.nepochsvec[obs];
  const double rotation = std::atan2(ev.p_y0[0], ev.p_x0[0]);
  double alpha = ev.alpha * TO_RAD;
  if (fixture.param.pllxMultiplyer)
    {
      parallax trajectory = ev.pllx[obs];
      trajectory.provide_observables_NE(ev.piEN, ev.piEE, ev.tE_r, ev.thE);
      alpha = std::atan2(trajectory.mulam_r, trajectory.mubet_r) - pi;
    }
  double tau = (ev.epoch[index] - ev.t0) / ev.tE_r;
  double impact = ev.u0;
  if (fixture.param.pllxMultiplyer)
    {
      tau -= ev.pllx[obs].Nshift[local] * ev.piEN +
        ev.pllx[obs].Eshift[local] * ev.piEE;
      impact += -ev.pllx[obs].Nshift[local] * ev.piEE +
        ev.pllx[obs].Eshift[local] * ev.piEN;
    }

  const double q = ev.p_q[0], s = ev.p_s0[0], rho = ev.rs;
  const double dx = s * std::cos(rotation), dy = s * std::sin(rotation);
  const double hostX = -q * dx / (1.0 + q);
  const double hostY = -q * dy / (1.0 + q);
  const double sourceX = impact * std::sin(alpha) + tau * std::cos(alpha) + hostX;
  const double sourceY = -impact * std::cos(alpha) + tau * std::sin(alpha) + hostY;
  double masses[2] = {1.0, q};
  VBcomplex lenses[2] = {VBcomplex(hostX, hostY), VBcomplex(hostX + dx, hostY + dy)};
  VBMicrolensing explicitSolver;
  explicitSolver.SetMethod(VBMicrolensing::Method::Multipoly);
  explicitSolver.Tol = explicitSolver.RelTol = 1e-10;
  explicitSolver.SetLensGeometry(2, masses, lenses);
  return explicitSolver.MultiMag(sourceX, sourceY, rho, 1e-10);
}

double explicitGeneratorMagnification(const event& ev, int index)
{
  double masses[2] = {1.0, ev.p_q[0]};
  VBcomplex lenses[2] = {
    VBcomplex(ev.xlens[0][index], ev.ylens[0][index]),
    VBcomplex(ev.xlens[1][index], ev.ylens[1][index])
  };
  VBMicrolensing explicitSolver;
  explicitSolver.SetMethod(VBMicrolensing::Method::Multipoly);
  explicitSolver.Tol = explicitSolver.RelTol = 1e-10;
  explicitSolver.SetLensGeometry(2, masses, lenses);
  return explicitSolver.MultiMag(ev.xsrc[0][index], ev.ysrc[0][index],
                                 ev.rs, 1e-10);
}

void checkGeneralGeneratorCOM()
{
  for (bool liveLensOrbit : {false, true})
    {
      Fixture fixture(0.05, true, liveLensOrbit);
      const event& ev = fixture.data;
      constexpr int index = 5;
      const double x = ev.xsrc[0][index], y = ev.ysrc[0][index];
      const double d0 = std::hypot(x - ev.xlens[0][index],
                                   y - ev.ylens[0][index]);
      const double d1 = std::hypot(x - ev.xlens[1][index],
                                   y - ev.ylens[1][index]);
      const double expected = explicitGeneratorMagnification(ev, index);
      require(std::min(d0, d1) > 10.0 * ev.rs,
              "general generator regression sample avoids both lens centers");
      require(close(ev.mu_src[0][index], expected, 2e-6, 2e-5),
              "general generator magnification matches explicit current lens geometry");
      std::cerr << "General generator COM orbit="
                << (liveLensOrbit ? "live" : "static")
                << " generated=" << ev.mu_src[0][index]
                << " MultiMag=" << expected << " dmin=" << std::min(d0, d1)
                << '\n';
    }
}

double oldPrimaryTransform(VBMicrolensing& solver, double q, double s,
                           double rotation, double trackX, double trackY,
                           double rho)
{
  const double a1 = -q * s / (1.0 + q);
  const double x = trackX + a1 * std::cos(rotation + pi);
  const double y = trackY + a1 * std::sin(rotation + pi);
  const double cr = std::cos(-rotation), sr = std::sin(-rotation);
  solver.a1 = 0.0;
  return solver.BinaryMag2(s, q, cr * x - sr * y, sr * x + cr * y, rho);
}

double explicitBinaryMagnification(double q, double s, double rotation,
                                   double trackX, double trackY, double rho)
{
  const double hostDistance = q * s / (1.0 + q);
  const double cr = std::cos(rotation), sr = std::sin(rotation);
  const double hostX = -hostDistance * cr, hostY = -hostDistance * sr;
  const double companionX = (s - hostDistance) * cr;
  const double companionY = (s - hostDistance) * sr;
  double masses[2] = {1.0, q};
  VBcomplex lenses[2] = {VBcomplex(hostX, hostY), VBcomplex(companionX, companionY)};
  VBMicrolensing explicitSolver;
  explicitSolver.SetMethod(VBMicrolensing::Method::Multipoly);
  explicitSolver.Tol = explicitSolver.RelTol = 1e-10;
  explicitSolver.SetLensGeometry(2, masses, lenses);
  return explicitSolver.MultiMag(hostX + trackX, hostY + trackY, rho, 1e-10);
}

double productionPrimaryTransform(VBMicrolensing& solver, double q, double s,
                                  double rotation, double trackX,
                                  double trackY, double rho)
{
  filekeywords param{};
  event ev{};
  ev.vbm = &solver;
  ev.p_q = {q};
  ev.p_s0 = {s};
  ev.p_x0 = {s * std::cos(rotation)};
  ev.p_y0 = {s * std::sin(rotation)};
  ev.tE_r = 1.0;
  ev.rs = rho;
  ev.u0 = -trackY;
  ev.epoch = {trackX};
  ev.obsidx = {0};
  DerivContext context{0, 0, &param, &ev, nullptr};
  return deriv_wrapper(0.0, &context);
}

void checkDenseBinaryCurves()
{
  // Compare the pre-fix primary-centered transform, the live production
  // deriv_wrapper, and explicit lens-plane Multipoly values over dense
  // finite-source tracks.
  struct Curve { double q, s, rotation, firstX, lastX; std::vector<double> y; };
  const std::vector<Curve> curves = {
    {0.001, 0.6, 0.0, -0.03, 0.03, {0.001, 0.01}},
    {0.1, 1.1, 1.1, -0.15, 0.15, {0.001, 0.01}},
    {1.0, 2.0, 0.0, -0.25, 0.25, {0.001, 0.01, 0.05}}
  };
  constexpr double step = 0.0005, rho = 0.001, tolerance = 1e-10;
  constexpr double strictAbs = 2e-6, strictRel = 2e-5;
  int total = 0;
  std::cerr << std::setprecision(8);
  for (const Curve& curve : curves)
    {
      VBMicrolensing oldSolver, productionSolver, polySolver;
      oldSolver.Tol = oldSolver.RelTol = tolerance;
      productionSolver.Tol = productionSolver.RelTol = tolerance;
      polySolver.SetMethod(VBMicrolensing::Method::Multipoly);
      polySolver.Tol = polySolver.RelTol = tolerance;
      const double cr = std::cos(curve.rotation), sr = std::sin(curve.rotation);
      const double primaryToCom = curve.q * curve.s / (1.0 + curve.q);
      double masses[2] = {1.0, curve.q};
      VBcomplex lenses[2] = {
        VBcomplex(-primaryToCom * cr, -primaryToCom * sr),
        VBcomplex((curve.s - primaryToCom) * cr,
                  (curve.s - primaryToCom) * sr)
      };
      polySolver.SetLensGeometry(2, masses, lenses);
      int samples = 0, allFinite = 0, oldNonfinite = 0;
      int productionNonfinite = 0, polyNonfinite = 0;
      int oldStrictMisses = 0, productionStrictMisses = 0, knownFiniteOutliers = 0;
      double oldMaxAbs = 0.0, oldMaxRel = 0.0;
      double productionMaxRel = 0.0, productionMaxAbs = 0.0;
      double maxX = 0.0, maxY = 0.0, maxCurrent = 0.0, maxPoly = 0.0;
      const int nSteps = static_cast<int>(std::round((curve.lastX - curve.firstX) / step));
      for (double y : curve.y)
        for (int i = 0; i <= nSteps; ++i)
          {
            const double x = curve.firstX + i * step;
            require(std::hypot(x, y) > 0.0 &&
                    std::hypot(x - curve.s * cr, y - curve.s * sr) > 0.0,
                    "dense trajectory does not pass through either lens center");
            const double oldValue = oldPrimaryTransform(oldSolver, curve.q, curve.s,
              curve.rotation, x, y, rho);
            const double currentValue = productionPrimaryTransform(productionSolver,
              curve.q, curve.s, curve.rotation, x, y, rho);
            const double polyValue = polySolver.MultiMag(lenses[0].re + x,
              lenses[0].im + y, rho, tolerance);
            ++samples;
            ++total;
            if (!std::isfinite(oldValue)) ++oldNonfinite;
            if (!std::isfinite(currentValue)) ++productionNonfinite;
            if (!std::isfinite(polyValue)) ++polyNonfinite;
            if (!std::isfinite(oldValue) || !std::isfinite(currentValue) ||
                !std::isfinite(polyValue))
              {
                std::cerr << "Dense solver nonfinite q=" << curve.q << " s="
                  << curve.s << " rot=" << curve.rotation << " x=" << x
                  << " y=" << y << " old=" << oldValue << " current="
                  << currentValue << " Multipoly=" << polyValue << '\n';
                continue;
              }
            ++allFinite;
            const double oldRelative = std::abs(oldValue - polyValue) /
              std::max(std::abs(oldValue), std::abs(polyValue));
            const double oldAbsolute = std::abs(oldValue - polyValue);
            const double currentAbsolute = std::abs(currentValue - polyValue);
            const double currentRelative = currentAbsolute /
              std::max(std::abs(currentValue), std::abs(polyValue));
            if (!close(oldValue, polyValue, strictAbs, strictRel)) ++oldStrictMisses;
            if (!close(currentValue, polyValue, strictAbs, strictRel))
              {
                ++productionStrictMisses;
                const bool knownCentralCausticOutlier = curve.q == 1.0 &&
                  curve.s == 2.0 && std::abs(x) < 1e-12 && y == 0.001;
                require(knownCentralCausticOutlier &&
                        close(currentValue, polyValue, strictAbs, 5e-4),
                        "only the recorded 0.047% central-caustic solver outlier exceeds strict tolerance");
                ++knownFiniteOutliers;
              }
            if (currentRelative > productionMaxRel)
              {
                productionMaxRel = currentRelative;
                productionMaxAbs = currentAbsolute;
                maxX = x; maxY = y; maxCurrent = currentValue; maxPoly = polyValue;
              }
            oldMaxAbs = std::max(oldMaxAbs, oldAbsolute);
            oldMaxRel = std::max(oldMaxRel, oldRelative);
          }
      // Established finite-source quadrature differences are measured and reported;
      // the strict misses remain visible instead of broadening the pass tolerance.
      std::cerr << "Dense curve q=" << curve.q << " s=" << curve.s
        << " rot=" << curve.rotation << " samples=" << samples
        << " all-finite=" << allFinite << " old-nonfinite=" << oldNonfinite
        << " current-nonfinite=" << productionNonfinite
        << " Multipoly-nonfinite=" << polyNonfinite
        << " old-strict-misses=" << oldStrictMisses
        << " current-strict-misses=" << productionStrictMisses
        << " known-finite-outliers=" << knownFiniteOutliers
        << " current-max-rel=" << productionMaxRel
        << " current-max-abs=" << productionMaxAbs
        << " current-max-at=(" << maxX << ',' << maxY << ") values=("
        << maxCurrent << ',' << maxPoly << ") old-max-rel=" << oldMaxRel
        << " old-max-abs=" << oldMaxAbs << '\n';
      require(samples == (curve.q == 0.001 ? 242 : curve.q == 0.1 ? 1202 : 3003),
              "dense finite-source curve has the requested sample count");
      require(oldStrictMisses > 0,
              "finite base-transform curve disagrees with explicit lenses");
      require(productionNonfinite == 0 && polyNonfinite == 0,
              "dense production and Multipoly curves remain finite");
    }
  require(total == 4447, "dense finite-source audit covers 4,447 samples");

  // This off-grid solver stress point is a known NaN in BinaryMag2. It is
  // deliberately classified separately, never used as evidence of the origin fix.
  VBMicrolensing oldSolver, productionSolver, polySolver;
  oldSolver.Tol = oldSolver.RelTol = productionSolver.Tol =
    productionSolver.RelTol = polySolver.Tol = polySolver.RelTol = tolerance;
  polySolver.SetMethod(VBMicrolensing::Method::Multipoly);
  double masses[2] = {1.0, 0.001};
  VBcomplex lenses[2] = {VBcomplex(-0.0006 / 1.001, 0.0),
                         VBcomplex(0.6 - 0.0006 / 1.001, 0.0)};
  polySolver.SetLensGeometry(2, masses, lenses);
  const double oldValue = oldPrimaryTransform(oldSolver, 0.001, 0.6, 0.0,
    0.00025, 0.001, rho);
  const double currentValue = productionPrimaryTransform(productionSolver,
    0.001, 0.6, 0.0, 0.00025, 0.001, rho);
  const double polyValue = polySolver.MultiMag(lenses[0].re + 0.00025,
    lenses[0].im + 0.001, rho, tolerance);
  std::cerr << "Known off-grid solver point q=0.001 s=0.6 rot=0 x=0.00025 y=0.001"
    << " old=" << oldValue << " current=" << currentValue
    << " Multipoly=" << polyValue << '\n';
  const bool expectedNan = std::isnan(currentValue) && std::isfinite(polyValue);
  const bool withinFiniteEnvelope = std::isfinite(currentValue) &&
    close(currentValue, polyValue, 2e-6, 5e-4);
  require(std::isfinite(oldValue) && std::isfinite(polyValue) &&
          (expectedNan || withinFiniteEnvelope),
          "known BinaryMag2 stress point is either nonfinite or within its recorded finite solver envelope");
}

void checkCenterOfMassSign()
{
  constexpr double q = 1.0, s = 2.0, trackX = 0.2, trackY = 0.31, rho = 0.001;
  for (double rotation : {0.0, 0.7})
    {
      VBMicrolensing oldSolver, productionSolver;
      oldSolver.Tol = oldSolver.RelTol = productionSolver.Tol = 1e-10;
      oldSolver.RelTol = productionSolver.RelTol = 1e-10;
      const double oldValue = oldPrimaryTransform(oldSolver, q, s, rotation,
                                                   trackX, trackY, rho);
      const double corrected = productionPrimaryTransform(productionSolver, q, s,
                                                           rotation, trackX, trackY, rho);
      const double reference = explicitBinaryMagnification(q, s, rotation,
                                                            trackX, trackY, rho);
      require(close(corrected, reference, 2e-6, 2e-5),
              "Fisher primary-relative track matches explicit binary lenses");
      require(std::abs(oldValue - reference) / reference > 0.1,
              "pre-fix CoM sign gives a finite, materially different magnification");
      std::cerr << "COM sign q=1 s=2 rot=" << rotation << " old=" << oldValue
                << " corrected=" << corrected << " Multipoly=" << reference << '\n';
    }
}

void checkParallaxStateAndDisabledParallax()
{
  // A piEN perturbation changes the stored observer shifts; the following t0
  // evaluation must reconstruct its own state rather than inherit that change.
  Fixture fixture;
  const int index = 16;
  DerivContext ordinary{0, index, &fixture.param, &fixture.data, fixture.world};
  DerivContext perturbed{7, index, &fixture.param, &fixture.data, fixture.world};
  const double before = deriv_wrapper(0.0, &ordinary);
  deriv_wrapper(0.04, &perturbed);
  const double after = deriv_wrapper(0.0, &ordinary);
  require(close(after, before, 1e-12, 1e-12),
          "a parameter perturbation does not leave parallax shifts for the next evaluation");

  Fixture noParallax(0.05, false);
  noParallax.data.alpha = 31.0;
  DerivContext noParallaxContext{0, index, &noParallax.param,
                                 &noParallax.data, noParallax.world};
  const double configured = deriv_wrapper(0.0, &noParallaxContext);
  noParallax.data.piEN = -0.8;
  noParallax.data.piEE = 0.6;
  const double unchanged = deriv_wrapper(0.0, &noParallaxContext);
  require(close(unchanged, configured, 1e-12, 1e-12),
          "disabled parallax uses the event trajectory angle and no parallax shifts");
}

void checkStaticFluxDerivatives()
{
  // The generated lightcurve includes live orbital motion, while Fisher flux
  // derivatives deliberately use the fixed projected lens configuration.
  Fixture fixture;
  event& ev = fixture.data;
  const std::vector<double> generatedMagnification = ev.Atrue;
  std::vector<double> staticMagnification(ev.nepochs);
  for (int i = 0; i < ev.nepochs; ++i)
    {
      DerivContext context{-1, i, &fixture.param, &ev, fixture.world};
      staticMagnification[i] = deriv_wrapper(0.0, &context);
      if (i % 11 == 0)
        require(close(staticMagnification[i], explicitLensMagnification(fixture, i),
                      2e-6, 2e-5),
                "static Fisher magnification matches explicit lenses");
      ev.Atrue[i] = 1.0 + ev.fs[ev.obsidx[i]] * (generatedMagnification[i] - 1.0);
    }

  fisherMatrix(&fixture.param, &ev, fixture.world, &fixture.sources, &fixture.lenses);
  bool liveAndStaticDiffer = false;
  for (int i = 0; i < ev.nepochs; ++i)
    {
      const int obs = ev.obsidx[i];
      const int base = 9 + 2 * obs;
      const double expectedBase = 1.0 + ev.fs[obs] * (staticMagnification[i] - 1.0);
      const double expectedFraction = staticMagnification[i] - 1.0;
      require(close(ev.dF[base * ev.nepochs + i], expectedBase, 1e-12, 1e-12),
              "Fbase derivative follows the fixed Fisher magnification");
      require(close(ev.dF[(base + 1) * ev.nepochs + i], expectedFraction,
                    1e-12, 1e-12),
              "fs derivative follows the fixed Fisher magnification");
      liveAndStaticDiffer |= std::abs(ev.Atrue[i] - expectedBase) > 1e-5;
    }
  require(liveAndStaticDiffer,
          "generated live-orbit flux differs from the fixed Fisher model at some epochs");
}

void checkInterpolationUsesComputedEpochs()
{
  // Force the existing suppression condition once and verify interpolation
  // reads computed same-observer neighbours rather than zero-initialized slots.
  Fixture fixture;
  event& ev = fixture.data;
  constexpr int parameter = 0;
  std::vector<double> raw(ev.nepochs);
  gsl_function function;
  function.function = deriv_wrapper;
  for (int i = 0; i < ev.nepochs; ++i)
    {
      DerivContext context{parameter, i, &fixture.param, &ev, fixture.world};
      function.params = &context;
      double result = 0.0, error = 0.0;
      gsl_deriv_central(&function, 0.0, 1e-4, &result, &error);
      raw[i] = result * ev.fs[ev.obsidx[i]];
    }

  const double expected = (raw[1] + raw[2] + raw[3]) / 3.0;
  require(std::isfinite(expected) && expected != 0.0,
          "valid simulated neighbors have a nonzero derivative");
  suppressFirstDerivative = true;
  fisherMatrix(&fixture.param, &ev, fixture.world, &fixture.sources, &fixture.lenses);
  suppressFirstDerivative = false;
  const double corrected = ev.dF[parameter * ev.nepochs];
  require(close(corrected, expected, 1e-7, 1e-6),
          "suppressed first epoch uses computed same-observer neighbors");
  require(corrected != 0.0,
          "old zero-initialized future neighbors would have produced zero");
  std::cerr << "First-epoch Fisher derivative: old=0, corrected=" << corrected << '\n';
}

} // namespace

int main()
{
  {
    QuietOutput quiet;
    checkGeneralGeneratorCOM();
    checkCenterOfMassSign();
    checkDenseBinaryCurves();
    checkParallaxStateAndDisabledParallax();
    checkStaticFluxDerivatives();
    checkInterpolationUsesComputedEpochs();
  }
  std::cout << "Fisher audit checks passed\n";
  return 0;
}
