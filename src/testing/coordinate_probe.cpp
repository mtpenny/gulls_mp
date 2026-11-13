#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "astroFns.h"
#include "coords.h"
#include "parallax.h"

int main()
{
  coords c;
  std::string operation;
  std::cout << std::setprecision(17);

  while(std::cin >> operation)
    {
      double x, y, z, a, b;
      std::cout << operation;

      if(operation == "ad2lb")
        {
          std::cin >> a >> b;
          std::vector<double> result;
          c.ad2lb(a, b, &x, &y);
          c.ad2lb(a, b, &result);
          std::cout << ' ' << x << ' ' << y << ' ' << result[0] << ' ' << result[1];
        }
      else if(operation == "lb2ad")
        {
          std::cin >> a >> b;
          std::vector<double> result;
          c.lb2ad(a, b, &x, &y);
          c.lb2ad(a, b, &result);
          std::cout << ' ' << x << ' ' << y << ' ' << result[0] << ' ' << result[1];
        }
      else if(operation == "ad2ecl")
        {
          std::cin >> a >> b;
          std::vector<double> result;
          c.ad2ecl(a, b, &x, &y);
          c.ad2ecl(a, b, &result);
          std::cout << ' ' << x << ' ' << y << ' ' << result[0] << ' ' << result[1];
        }
      else if(operation == "ecl2ad")
        {
          std::cin >> a >> b;
          std::vector<double> result;
          c.ecl2ad(a, b, &x, &y);
          c.ecl2ad(a, b, &result);
          std::cout << ' ' << x << ' ' << y << ' ' << result[0] << ' ' << result[1];
        }
      else if(operation == "muad2lb")
        {
          std::cin >> a >> b >> x >> y;
          double first, second;
          c.muad2lb(a, b, x, y, &first, &second);
          std::cout << ' ' << first << ' ' << second;
        }
      else if(operation == "mulb2ad")
        {
          std::cin >> a >> b >> x >> y;
          double first, second;
          c.mulb2ad(a, b, x, y, &first, &second);
          std::cout << ' ' << first << ' ' << second;
        }
      else if(operation == "muad2ecl")
        {
          std::cin >> a >> b >> x >> y;
          double first, second;
          c.muad2ecl(a, b, x, y, &first, &second);
          std::cout << ' ' << first << ' ' << second;
        }
      else if(operation == "muecl2ad")
        {
          std::cin >> a >> b >> x >> y;
          double first, second;
          c.muecl2ad(a, b, x, y, &first, &second);
          std::cout << ' ' << first << ' ' << second;
        }
      else if(operation == "xyz2rlb")
        {
          std::cin >> x >> y >> z;
          const double input_x = x, input_y = y, input_z = z;
          std::vector<double> result;
          c.XYZ2Rlb(input_x, input_y, input_z, &a, &b, &z);
          c.XYZ2Rlb(std::vector<double>{input_x, input_y, input_z}, &result);
          std::cout << ' ' << a << ' ' << b << ' ' << z
                    << ' ' << result[0] << ' ' << result[1] << ' ' << result[2];
        }
      else if(operation == "rlb2xyz")
        {
          std::cin >> x >> y >> z;
          const double input_r = x, input_l = y, input_b = z;
          std::vector<double> result;
          c.Rlb2XYZ(input_r, input_l, input_b, &x, &y, &z);
          c.Rlb2XYZ(std::vector<double>{input_r, input_l, input_b}, &result);
          std::cout << ' ' << x << ' ' << y << ' ' << z
                    << ' ' << result[0] << ' ' << result[1] << ' ' << result[2];
        }
      else if(operation == "uvw")
        {
          std::cin >> a >> b >> x >> y >> z;
          double w;
          std::cin >> w;
          std::vector<double> result;
          c.uvw2rvmuadlb(a, b, x, y, z, w, &result);
          for(double value : result) std::cout << ' ' << value;
        }
      else if(operation == "eq2gal")
        {
          char direction;
          std::cin >> direction >> a >> b;
          eq2gal(a, b, direction, &x, &y);
          std::cout << ' ' << x << ' ' << y;
        }
      else if(operation == "set_lb")
        {
          std::cin >> a >> b;
          parallax p;
          p.set_lb(a, b, 1);
          std::cout << ' ' << p.a << ' ' << p.d << ' ' << p.lam << ' ' << p.bet;
        }
      else if(operation == "lst")
        {
          double jd, longitude;
          std::cin >> jd >> longitude;
          lst(jd, longitude, &x);
          std::cout << ' ' << x;
        }
      else if(operation == "geod2geoc")
        {
          double longitude, latitude, height;
          std::string ellipsoid;
          std::cin >> longitude >> latitude >> height >> ellipsoid;
          double geocentric[3], cartesian[3];
          geod2geoc(longitude, latitude, height, ellipsoid.c_str(), geocentric, cartesian);
          for(double value : geocentric) std::cout << ' ' << value;
          for(double value : cartesian) std::cout << ' ' << value;
        }
      else if(operation == "eq2horiz")
        {
          char direction;
          double jd, latitude, longitude;
          std::cin >> direction >> a >> b >> jd >> latitude >> longitude;
          eq2horiz(a, b, jd, latitude, longitude, direction, &x, &y);
          std::cout << ' ' << x << ' ' << y;
        }
      else if(operation == "eq2eclip")
        {
          std::cin >> a >> b;
          eq2eclip(a, b, &x, &y);
          std::cout << ' ' << x << ' ' << y;
        }
      else if(operation == "gal2eclip")
        {
          std::cin >> a >> b;
          gal2eclip(a, b, &x, &y);
          std::cout << ' ' << x << ' ' << y;
        }
      else if(operation == "units")
        {
          std::cin >> x;
          std::cout << ' ' << c.rad2hrs(x) << ' ' << c.hrs2rad(x)
                    << ' ' << c.rad2deg(x) << ' ' << c.deg2rad(x);
        }
      else if(operation == "fold")
        {
          std::cin >> x >> y >> z;
          std::cout << ' ' << c.fold(x, y, z);
        }
      else
        {
          std::cerr << "unknown operation: " << operation << '\n';
          return 2;
        }
      std::cout << '\n';
    }
}
