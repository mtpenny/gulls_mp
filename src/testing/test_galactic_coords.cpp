#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "astroFns.h"
#include "parallax.h"

namespace
{
  constexpr double test_pi = 3.14159265358979323846;
  constexpr double mas_per_rad = 180.0 / test_pi * 3600.0 * 1000.0;

  double wrap(double angle)
  {
    return std::remainder(angle, 2.0*test_pi);
  }

  double separation_mas(double lon1, double lat1, double lon2, double lat2)
  {
    const double c1 = std::cos(lat1), c2 = std::cos(lat2);
    const double x1 = c1*std::cos(lon1), y1 = c1*std::sin(lon1), z1 = std::sin(lat1);
    const double x2 = c2*std::cos(lon2), y2 = c2*std::sin(lon2), z2 = std::sin(lat2);
    const double cx = y1*z2 - z1*y2, cy = z1*x2 - x1*z2, cz = x1*y2 - y1*x2;
    const double dot = x1*x2 + y1*y2 + z1*z2;
    return std::atan2(std::sqrt(cx*cx + cy*cy + cz*cz), dot)*mas_per_rad;
  }

  void require(bool condition, const std::string& message)
  {
    if(!condition) throw std::runtime_error(message);
  }

  struct GalacticReference
  {
    double l_deg, b_deg, ra_fk5_deg, dec_fk5_deg;
    double mu_l, mu_b, mu_ra_fk5, mu_dec_fk5;
  };

  // Astropy 7.2.0 FK5(equinox=J2000) reference positions and proper motions.
  // The first Galactic direction is inferred from the reported table RA/Dec.
  const GalacticReference references[] = {
    {0.3224581385681185, -1.8869917206356943,
     268.451164329658, -29.628745429084, 3.2, -4.1, 5.155682639424, 0.684789399408},
    {359.8, -0.1, 266.383408768680, -29.158981297122,
     -2.3, 5.7, -6.063357050489, 1.007820062456},
    {123.0, 70.0, 192.825256350203, 47.128233084048,
     1.25, -0.75, -1.251164259625, 0.748056144576},
    {203.123, 89.9, 192.748745868932, 27.145244041119,
     6.5, 2.2, 1.066456503552, -6.778839910046}
  };

  void check_galactic_derivative(coords& c, double l, double b, double mul, double mub,
                                  double mua, double mud)
  {
    const double dt = 0.05;
    double a, d, ap, dp, am, dm;
    c.lb2ad(l, b, &a, &d);
    const double dl = mul*dt/(mas_per_rad*std::cos(b));
    const double db = mub*dt/mas_per_rad;
    c.lb2ad(l + dl, b + db, &ap, &dp);
    c.lb2ad(l - dl, b - db, &am, &dm);
    const double fd_ra = wrap(ap-am)*std::cos(d)*mas_per_rad/(2.0*dt);
    const double fd_dec = (dp-dm)*mas_per_rad/(2.0*dt);
    require(std::abs(fd_ra-mua) < 2e-4 && std::abs(fd_dec-mud) < 2e-4,
            "Galactic-to-equatorial proper motion disagrees with position derivative");
  }

  void check_equatorial_derivative(coords& c, double a, double d, double mua, double mud,
                                    double mul, double mub)
  {
    const double dt = 0.05;
    double l, b, lp, bp, lm, bm;
    c.ad2lb(a, d, &l, &b);
    const double da = mua*dt/(mas_per_rad*std::cos(d));
    const double dd = mud*dt/mas_per_rad;
    c.ad2lb(a + da, d + dd, &lp, &bp);
    c.ad2lb(a - da, d - dd, &lm, &bm);
    const double fd_l = wrap(lp-lm)*std::cos(b)*mas_per_rad/(2.0*dt);
    const double fd_b = (bp-bm)*mas_per_rad/(2.0*dt);
    require(std::abs(fd_l-mul) < 2e-4 && std::abs(fd_b-mub) < 2e-4,
            "Equatorial-to-Galactic proper motion disagrees with position derivative");
  }
}

int main()
{
  coords c;
  double max_fk5_separation_mas = 0.0;
  double max_pm_error = 0.0;
  double event_path_separation_mas = 0.0;

  for(const auto& ref : references)
    {
      const double l = ref.l_deg*test_pi/180.0, b = ref.b_deg*test_pi/180.0;
      double a, d, table_a, table_d, l_back, b_back;
      c.lb2ad(l, b, &a, &d);
      eq2gal(l, b, 'g', &table_a, &table_d);
      if(table_a < 0.0) table_a += 2.0*test_pi; // Match the production caller's RA wrap.
      const double path_sep = separation_mas(a, d, table_a, table_d);
      require(path_sep < 1e-5, "table and parallax Galactic conversions diverged");
      if(ref.l_deg < 1.0) event_path_separation_mas = path_sep;

      const double fk5_a = ref.ra_fk5_deg*test_pi/180.0, fk5_d = ref.dec_fk5_deg*test_pi/180.0;
      const double fk5_sep = separation_mas(a, d, fk5_a, fk5_d);
      max_fk5_separation_mas = std::max(max_fk5_separation_mas, fk5_sep);
      // Allow for the published matrix's 10-decimal coefficients and frame realization.
      require(fk5_sep < 20.0, "position differs from independent FK5/J2000 reference by over 20 mas");

      eq2gal(a, d, 'e', &l_back, &b_back);
      if(l_back < 0.0) l_back += 2.0*test_pi;
      const double inverse_sep = separation_mas(l, b, l_back, b_back);
      // The inverse uses the transpose of the rounded rotation matrix.
      require(inverse_sep < 0.1,
              "equatorial/Galactic inverse conversion failed");

      double mua, mud, mul_back, mub_back;
      c.mulb2ad(l, b, ref.mu_l, ref.mu_b, &mua, &mud);
      const double pm_error = std::max(std::abs(mua-ref.mu_ra_fk5), std::abs(mud-ref.mu_dec_fk5));
      max_pm_error = std::max(max_pm_error, pm_error);
      require(pm_error < 2e-5, "proper motion differs from independent FK5/J2000 reference");
      require(std::abs(std::hypot(mua, mud) - std::hypot(ref.mu_l, ref.mu_b)) < 1e-6,
              "proper-motion rotation did not preserve vector norm");
      c.muad2lb(a, d, mua, mud, &mul_back, &mub_back);
      require(std::abs(mul_back-ref.mu_l) < 1e-7 && std::abs(mub_back-ref.mu_b) < 1e-7,
              "proper-motion forward/reverse transform failed");
      check_galactic_derivative(c, l, b, ref.mu_l, ref.mu_b, mua, mud);
    }

  parallax pllx;
  const auto& event = references[0];
  pllx.set_lb(event.l_deg, event.b_deg);
  double event_table_a, event_table_d;
  eq2gal(event.l_deg*test_pi/180.0, event.b_deg*test_pi/180.0, 'g', &event_table_a, &event_table_d);
  if(event_table_a < 0.0) event_table_a += 2.0*test_pi;
  require(separation_mas(pllx.a, pllx.d, event_table_a, event_table_d) < 1e-5,
          "parallax::set_lb does not use the published Galactic conversion");

  const double equatorial_cases[][4] = {
    {359.999, 12.3, 5.2, -1.1}, {0.001, -45.0, -2.3, 4.4},
    {42.0, 89.8, 3.1, -2.7}, {219.0, -89.8, -1.6, 3.8}
  };
  for(const auto& test : equatorial_cases)
    {
      const double a = test[0]*test_pi/180.0, d = test[1]*test_pi/180.0;
      double l, b, a_back, d_back, mul, mub, mua_back, mud_back;
      c.ad2lb(a, d, &l, &b);
      c.lb2ad(l, b, &a_back, &d_back);
      require(separation_mas(a, d, a_back, d_back) < 0.1,
              "position inverse failed near RA wrap or a celestial pole");
      c.muad2lb(a, d, test[2], test[3], &mul, &mub);
      require(std::abs(std::hypot(mul, mub) - std::hypot(test[2], test[3])) < 1e-6,
              "reverse proper-motion rotation did not preserve vector norm");
      c.mulb2ad(l, b, mul, mub, &mua_back, &mud_back);
      require(std::abs(mua_back-test[2]) < 1e-7 && std::abs(mud_back-test[3]) < 1e-7,
              "proper-motion inverse failed near RA wrap or a celestial pole");
      check_equatorial_derivative(c, a, d, test[2], test[3], mul, mub);
    }

  std::cout << "event table/parallax separation: " << event_path_separation_mas/1000.0
            << " arcsec; max FK5/J2000 position offset: " << max_fk5_separation_mas/1000.0
            << " arcsec; max proper-motion reference error: " << max_pm_error << " mas/yr\n";
}
