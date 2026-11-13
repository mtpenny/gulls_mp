#ifndef GALACTIC_ROTATION_H
#define GALACTIC_ROTATION_H

namespace galactic_rotation
{
  // Keep published J2000 RA/Dec and parallax helpers on one Galactic rotation.
  // Do not substitute a different equatorial realization in only one path.
  constexpr double equatorial_to_galactic[3][3] = {
    {-0.0548755604, -0.8734370902, -0.4838350155},
    { 0.4941094279, -0.4448296300,  0.7469822445},
    {-0.8676661490, -0.1980763734,  0.4559837762}
  };
}

#endif
