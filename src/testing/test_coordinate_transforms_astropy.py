"""Compare the production C++ coordinate APIs with explicit Astropy frames.

Run from the repository root with:
  python3 -m pip install -r src/testing/requirements_coordinates.txt
  python3 src/testing/test_coordinate_transforms_astropy.py
"""

import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

import astropy.units as u
import numpy as np
from astropy.coordinates import (
    AltAz,
    BarycentricMeanEcliptic,
    CartesianDifferential,
    CartesianRepresentation,
    EarthLocation,
    FK5,
    Galactic,
    HADec,
    SkyCoord,
)
from astropy.time import Time
from astropy.utils import iers


ROOT = Path(__file__).resolve().parents[2]
PI = math.pi
POSITION_FK5_TOL = 10 * u.mas
ECLIPTIC_FRAME_TOL = 35 * u.mas
PM_TOL = 2e-5 * u.mas / u.yr

# First l,b is reconstructed from the supplied published table RA/Dec; it is
# a regression direction, not a claim that these are the historical event inputs.
GALACTIC_DIRECTIONS_DEG = [
    (0.3224581385681185, -1.8869917206356943),
    (359.8, -0.1),
    (123.0, 70.0),
    (203.123, 89.9),
    (33.0, -89.9),
    (0.0, 0.0),
    (90.0, 0.0),
    (180.0, -30.0),
    (359.999, 0.1),
    (122.93191999727053, 27.128250001902302),
]


def _radians(degrees):
    return math.radians(degrees)


def _normalized_ra(radians):
    return radians % (2.0 * PI)


def _coordinate_separation(first, second):
    return first.separation(second)


def _probe_commands(commands):
    return [" ".join(str(value) for value in command) for command in commands]


class CoordinateTransformAstropyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temporary = tempfile.TemporaryDirectory(
            prefix=".coordinate-transform-test-", dir=ROOT
        )
        cls.probe_path = Path(cls._temporary.name) / "coordinate_probe"
        compiler = shlex.split(os.environ.get("CXX", "c++"))
        command = compiler + ["-std=c++17"]
        command += [
            f"-I{ROOT / 'src'}",
            f"-I{ROOT / 'src/headers'}",
            f"-I{ROOT / 'src/classes'}",
        ]

        # Homebrew keeps GSL keg-only; Ubuntu CI uses standard include/library paths.
        if sys.platform == "darwin" and shutil.which("brew"):
            prefix = subprocess.run(
                ["brew", "--prefix", "gsl"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
            command += [f"-I{prefix}/include", f"-L{prefix}/lib"]

        command += [
            str(ROOT / "src/testing/coordinate_probe.cpp"),
            str(ROOT / "src/astroFns.cpp"),
            str(ROOT / "src/classes/coords.cpp"),
            str(ROOT / "src/classes/parallax.cpp"),
            str(ROOT / "src/classes/ephem.cpp"),
            "-lgsl",
            "-lgslcblas",
            "-o",
            str(cls.probe_path),
        ]
        subprocess.run(command, cwd=ROOT, check=True)
        iers.conf.auto_download = False
        cls.fk5 = FK5(equinox=Time("J2000"))
        cls.ecliptic = BarycentricMeanEcliptic(equinox=Time("J2000"))

    @classmethod
    def tearDownClass(cls):
        cls._temporary.cleanup()

    def _probe(self, commands):
        lines = _probe_commands(commands)
        result = subprocess.run(
            [str(type(self).probe_path)],
            input="\n".join(lines) + "\n",
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        rows = [line.split() for line in result.stdout.splitlines()]
        self.assertEqual([row[0] for row in rows], [line.split()[0] for line in lines])
        return [np.asarray([float(value) for value in row[1:]]) for row in rows]

    def test_galactic_equatorial_positions_and_published_table_path(self):
        directions = GALACTIC_DIRECTIONS_DEG
        lb_commands = [["lb2ad", _radians(l), _radians(b)] for l, b in directions]
        lb_results = self._probe(lb_commands)
        eq_commands = []
        for result in lb_results:
            eq_commands.append(["ad2lb", result[0] % (2 * PI), result[1]])
        eq_results = self._probe(eq_commands)

        table_commands = [
            ["eq2gal", "g", _radians(l), _radians(b)] for l, b in directions
        ]
        table_results = self._probe(table_commands)
        table_inverse_results = self._probe(
            [["eq2gal", "e", _normalized_ra(result[0]), result[1]] for result in lb_results]
        )
        for (l, b), lb_result, eq_result, table_result, table_inverse in zip(
            directions, lb_results, eq_results, table_results, table_inverse_results
        ):
            # Both pointer and vector overloads exercise the same production path.
            np.testing.assert_allclose(lb_result[:2], lb_result[2:4], atol=2e-15)
            np.testing.assert_allclose(eq_result[:2], eq_result[2:4], atol=2e-15)

            expected_gal = SkyCoord(l=l * u.deg, b=b * u.deg, frame=Galactic())
            expected_fk5 = expected_gal.transform_to(self.fk5)
            actual_fk5 = SkyCoord(
                ra=_normalized_ra(lb_result[0]) * u.rad,
                dec=lb_result[1] * u.rad,
                frame=self.fk5,
            )
            self.assertLess(
                _coordinate_separation(actual_fk5, expected_fk5),
                POSITION_FK5_TOL,
            )

            actual_gal = SkyCoord(
                l=eq_result[0] * u.rad,
                b=eq_result[1] * u.rad,
                frame=Galactic(),
            )
            self.assertLess(
                _coordinate_separation(actual_gal, expected_gal),
                POSITION_FK5_TOL,
            )

            # The existing eq2gal table matrix remains an independent production path.
            actual_table = SkyCoord(
                ra=_normalized_ra(table_result[0]) * u.rad,
                dec=table_result[1] * u.rad,
                frame=self.fk5,
            )
            self.assertLess(
                _coordinate_separation(actual_table, expected_fk5),
                POSITION_FK5_TOL,
            )
            self.assertLess(
                _coordinate_separation(actual_table, actual_fk5),
                0.02 * u.mas,
            )
            inverse_table_gal = SkyCoord(
                l=table_inverse[0] * u.rad,
                b=table_inverse[1] * u.rad,
                frame=Galactic(),
            )
            self.assertLess(
                _coordinate_separation(inverse_table_gal, expected_gal),
                0.1 * u.mas,
            )

            # Round trips remain accurate through the rounded published matrix,
            # including both Galactic poles and the equatorial north-pole direction.
            self.assertLess(
                _coordinate_separation(actual_gal, expected_gal),
                POSITION_FK5_TOL,
            )

        event_l, event_b = directions[0]
        table = table_results[0]
        saved_table = SkyCoord(
            ra=268.4511629368 * u.deg,
            dec=-29.62874695772 * u.deg,
            frame=self.fk5,
        )
        actual_saved = SkyCoord(
            ra=_normalized_ra(table[0]) * u.rad,
            dec=table[1] * u.rad,
            frame=self.fk5,
        )
        self.assertLess(
            _coordinate_separation(actual_saved, saved_table),
            0.02 * u.mas,
        )

        # eq2gal's inverse direction accepts the FK5/J2000 table coordinates.
        inverse = self._probe(
            [["eq2gal", "e", _normalized_ra(table[0]), table[1]]]
        )[0]
        inverse_gal = SkyCoord(l=inverse[0] * u.rad, b=inverse[1] * u.rad, frame=Galactic())
        expected_gal = SkyCoord(
            l=event_l * u.deg, b=event_b * u.deg, frame=Galactic()
        )
        self.assertLess(_coordinate_separation(inverse_gal, expected_gal), 0.02 * u.mas)

        set_lb = self._probe([["set_lb", _radians(event_l), _radians(event_b)]])[0]
        self.assertLess(
            _coordinate_separation(
                SkyCoord(ra=_normalized_ra(set_lb[0]) * u.rad, dec=set_lb[1] * u.rad, frame=self.fk5),
                actual_saved,
            ),
            0.02 * u.mas,
        )
        expected_ecliptic = self._probe([["ad2ecl", set_lb[0], set_lb[1]]])[0]
        self.assertLess(
            _coordinate_separation(
                SkyCoord(lon=set_lb[2] * u.rad, lat=set_lb[3] * u.rad, frame=self.ecliptic),
                SkyCoord(lon=expected_ecliptic[0] * u.rad,
                         lat=expected_ecliptic[1] * u.rad,
                         frame=self.ecliptic),
            ),
            0.02 * u.mas,
        )

    def test_equatorial_ecliptic_positions_and_astrofn_helpers(self):
        directions = GALACTIC_DIRECTIONS_DEG
        gal_positions = self._probe(
            [["lb2ad", _radians(l), _radians(b)] for l, b in directions]
        )
        for (l, b), position in zip(directions, gal_positions):
            ra = _normalized_ra(position[0])
            fk_position = SkyCoord(ra=ra * u.rad, dec=position[1] * u.rad, frame=self.fk5)
            expected = fk_position.transform_to(self.ecliptic)
            coords_result = self._probe([["ad2ecl", ra, position[1]]])[0]
            np.testing.assert_allclose(coords_result[:2], coords_result[2:4], atol=2e-15)
            actual = SkyCoord(
                lon=coords_result[0] * u.rad,
                lat=coords_result[1] * u.rad,
                frame=self.ecliptic,
            )
            self.assertLess(_coordinate_separation(actual, expected), ECLIPTIC_FRAME_TOL)

            astrofn_result = self._probe(
                [["eq2eclip", math.degrees(ra), math.degrees(position[1])]]
            )[0]
            astrofn_coord = SkyCoord(
                lon=astrofn_result[0] * u.deg,
                lat=astrofn_result[1] * u.deg,
                frame=self.ecliptic,
            )
            self.assertLess(_coordinate_separation(astrofn_coord, actual), 0.02 * u.mas)

            galactic_ecliptic = self._probe([["gal2eclip", l, b]])[0]
            galactic_ecliptic_coord = SkyCoord(
                lon=galactic_ecliptic[0] * u.deg,
                lat=galactic_ecliptic[1] * u.deg,
                frame=self.ecliptic,
            )
            self.assertLess(
                _coordinate_separation(galactic_ecliptic_coord, expected),
                ECLIPTIC_FRAME_TOL,
            )
            self.assertLess(
                _coordinate_separation(galactic_ecliptic_coord, actual),
                0.05 * u.mas,
            )

        ecliptic_directions = [
            (10.0, 20.0),
            (359.999, 0.1),
            (123.0, 70.0),
            (45.0, -80.0),
            (10.0, 89.9),
            (230.0, -89.9),
        ]
        for longitude, latitude in ecliptic_directions:
            expected_ecliptic = SkyCoord(
                lon=longitude * u.deg,
                lat=latitude * u.deg,
                frame=self.ecliptic,
            )
            expected_fk5 = expected_ecliptic.transform_to(self.fk5)
            result = self._probe(
                [["ecl2ad", _radians(longitude), _radians(latitude)]]
            )[0]
            np.testing.assert_allclose(result[:2], result[2:4], atol=2e-15)
            actual_fk5 = SkyCoord(
                ra=_normalized_ra(result[0]) * u.rad,
                dec=result[1] * u.rad,
                frame=self.fk5,
            )
            self.assertLess(
                _coordinate_separation(actual_fk5, expected_fk5),
                ECLIPTIC_FRAME_TOL,
            )

    def test_galactic_and_ecliptic_proper_motion_rotations(self):
        galactic_pm = [
            (3.2, -4.1),
            (-2.3, 5.7),
            (6.5, 2.2),
            (1.25, -0.75),
            (-1.6, 3.8),
        ]
        gal_positions = self._probe(
            [["lb2ad", _radians(l), _radians(b)] for l, b in GALACTIC_DIRECTIONS_DEG]
        )
        for (l, b), position, (mu_l, mu_b) in zip(
            GALACTIC_DIRECTIONS_DEG, gal_positions, galactic_pm
        ):
            actual_from_gal = self._probe(
                [["mulb2ad", _radians(l), _radians(b), mu_l, mu_b]]
            )[0]
            expected_gal_star = SkyCoord(
                l=l * u.deg,
                b=b * u.deg,
                distance=1 * u.kpc,
                pm_l_cosb=mu_l * u.mas / u.yr,
                pm_b=mu_b * u.mas / u.yr,
                frame=Galactic(),
            )
            expected_fk5_star = expected_gal_star.transform_to(self.fk5)
            np.testing.assert_allclose(
                actual_from_gal,
                [expected_fk5_star.pm_ra_cosdec.value, expected_fk5_star.pm_dec.value],
                atol=PM_TOL.to_value(u.mas / u.yr),
            )

            mu_ra, mu_dec = 2.7 - mu_l / 3.0, -1.4 + mu_b / 4.0
            fk5_star = SkyCoord(
                ra=_normalized_ra(position[0]) * u.rad,
                dec=position[1] * u.rad,
                distance=1 * u.kpc,
                pm_ra_cosdec=mu_ra * u.mas / u.yr,
                pm_dec=mu_dec * u.mas / u.yr,
                frame=self.fk5,
            )
            if abs(b) < 89.0:
                # At |b|=89.9 deg, the 5-6 mas frame-realization offset changes
                # the longitude tangent basis by enough to shift components by
                # 3.4e-5 mas/yr. Near-pole behavior is instead checked by C++
                # position derivatives, PM norm preservation, and inverse tests.
                actual_fk5 = self._probe(
                    [["muad2lb", _normalized_ra(position[0]), position[1], mu_ra, mu_dec]]
                )[0]
                expected_gal_star = fk5_star.transform_to(Galactic())
                np.testing.assert_allclose(
                    actual_fk5,
                    [expected_gal_star.pm_l_cosb.value, expected_gal_star.pm_b.value],
                    atol=PM_TOL.to_value(u.mas / u.yr),
                )
            self.assertAlmostEqual(
                math.hypot(*actual_from_gal), math.hypot(mu_l, mu_b), delta=1e-9
            )

            actual_eclip_pm = self._probe(
                [["muad2ecl", _normalized_ra(position[0]), position[1], mu_ra, mu_dec]]
            )[0]
            expected_eclip_star = fk5_star.transform_to(self.ecliptic)
            np.testing.assert_allclose(
                actual_eclip_pm,
                [expected_eclip_star.pm_lon_coslat.value, expected_eclip_star.pm_lat.value],
                atol=PM_TOL.to_value(u.mas / u.yr),
            )

        for longitude, latitude in [
            (10.0, 20.0),
            (359.999, 0.1),
            (123.0, 70.0),
            (45.0, -80.0),
            (10.0, 89.9),
            (230.0, -89.9),
        ]:
            equatorial = self._probe(
                [["ecl2ad", _radians(longitude), _radians(latitude)]]
            )[0]
            ra, dec = _normalized_ra(equatorial[0]), equatorial[1]
            input_pm = (3.2, -4.1)
            fk5_star = SkyCoord(
                ra=ra * u.rad,
                dec=dec * u.rad,
                distance=1 * u.kpc,
                pm_ra_cosdec=input_pm[0] * u.mas / u.yr,
                pm_dec=input_pm[1] * u.mas / u.yr,
                frame=self.fk5,
            )
            expected_eclip = fk5_star.transform_to(self.ecliptic)
            actual = self._probe(
                [[
                    "muecl2ad",
                    ra,
                    dec,
                    expected_eclip.pm_lon_coslat.value,
                    expected_eclip.pm_lat.value,
                ]]
            )[0]
            if abs(latitude) < 89.0:
                np.testing.assert_allclose(actual, input_pm, atol=1e-5)
            else:
                # FK5/ICRS realization offsets rotate the tangent basis near the
                # ecliptic pole; retain the direct C++ inverse/norm check there.
                cpp_eclip_pm = self._probe(
                    [["muad2ecl", ra, dec, input_pm[0], input_pm[1]]]
                )[0]
                inverse = self._probe(
                    [["muecl2ad", ra, dec, cpp_eclip_pm[0], cpp_eclip_pm[1]]]
                )[0]
                np.testing.assert_allclose(inverse, input_pm, atol=1e-7)
                self.assertAlmostEqual(
                    math.hypot(*cpp_eclip_pm), math.hypot(*input_pm), delta=1e-9
                )

    def test_xyz_spherical_overloads_against_astropy(self):
        xyz_values = [(2.0, -1.0, 3.0), (-4.0, 2.0, -0.5), (0.0, -2.0, 1.5)]
        for x, y, z in xyz_values:
            actual = self._probe([["xyz2rlb", x, y, z]])[0]
            np.testing.assert_allclose(actual[:3], actual[3:6], atol=2e-15)
            expected = SkyCoord(
                CartesianRepresentation(x=x * u.kpc, y=y * u.kpc, z=z * u.kpc),
                frame=Galactic(),
            )
            self.assertAlmostEqual(actual[0], expected.distance.to_value(u.kpc), delta=1e-12)
            actual_spherical = SkyCoord(
                l=actual[1] * u.rad,
                b=actual[2] * u.rad,
                distance=actual[0] * u.kpc,
                frame=Galactic(),
            )
            self.assertLess(_coordinate_separation(actual_spherical, expected), 1e-7 * u.mas)

        spherical_values = [
            (4.0, _radians(359.8), _radians(-0.1)),
            (2.5, _radians(123.0), _radians(70.0)),
            (8.0, _radians(33.0), _radians(-89.9)),
        ]
        for radius, longitude, latitude in spherical_values:
            actual = self._probe([["rlb2xyz", radius, longitude, latitude]])[0]
            np.testing.assert_allclose(actual[:3], actual[3:6], atol=2e-15)
            expected = SkyCoord(
                l=longitude * u.rad,
                b=latitude * u.rad,
                distance=radius * u.kpc,
                frame=Galactic(),
            ).cartesian
            np.testing.assert_allclose(
                actual[:3],
                [expected.x.to_value(u.kpc), expected.y.to_value(u.kpc), expected.z.to_value(u.kpc)],
                atol=2e-14,
            )

    def test_uvw_velocity_projection_against_astropy(self):
        cases = [
            (12.0, 20.0, 4.0, 10.0, 20.0, 30.0),
            (359.8, -0.1, 1.2, -120.0, 30.0, 45.0),
            (220.0, 70.0, 8.0, 15.0, -80.0, 12.0),
            (203.123, 89.9, 2.0, 1.0, 2.0, 3.0),
        ]
        astropy_mas_per_year_per_kms_kpc = (
            (1 * u.km / u.s / u.kpc)
            .to(u.mas / u.yr, equivalencies=u.dimensionless_angles())
            .value
        )
        # coords.cpp retains its documented historical 4.74057 conversion.
        legacy_scale = (1.0 / 4.74057) / astropy_mas_per_year_per_kms_kpc
        for l, b, distance, U, V, W in cases:
            actual = self._probe(
                [["uvw", _radians(l), _radians(b), distance, U, V, W]]
            )[0]
            position = SkyCoord(
                l=l * u.deg,
                b=b * u.deg,
                distance=distance * u.kpc,
                frame=Galactic(),
            )
            differential = CartesianDifferential(
                d_x=U * u.km / u.s,
                d_y=V * u.km / u.s,
                d_z=W * u.km / u.s,
            )
            moving = SkyCoord(
                position.cartesian.with_differentials(differential), frame=Galactic()
            )
            expected_fk5 = moving.transform_to(self.fk5)
            expected_galactic = moving
            self.assertAlmostEqual(
                actual[0], expected_fk5.radial_velocity.to_value(u.km / u.s), delta=1e-6
            )
            np.testing.assert_allclose(
                actual[1:3],
                [
                    expected_fk5.pm_ra_cosdec.value * legacy_scale,
                    expected_fk5.pm_dec.value * legacy_scale,
                ],
                atol=5e-5,
            )
            np.testing.assert_allclose(
                actual[3:5],
                [
                    expected_galactic.pm_l_cosb.value * legacy_scale,
                    expected_galactic.pm_b.value * legacy_scale,
                ],
                atol=5e-5,
            )

    def test_astrofn_horizontal_sidereal_and_geodetic_helpers(self):
        jd = 2451545.0
        horizontal_cases = [
            (180.0, 20.0, 0.5, -1.2),
            (270.0, -30.0, 0.7, 2.1),
            (15.0, 45.0, -0.4, 0.8),
            (100.0, 5.0, 0.1, -2.2),
        ]
        for ra_deg, dec_deg, latitude, longitude in horizontal_cases:
            ra, dec = _radians(ra_deg), _radians(dec_deg)
            lst = self._probe([["lst", jd, longitude]])[0][0]
            time = Time(jd, format="jd", scale="ut1")
            location = EarthLocation.from_geodetic(
                lon=longitude * u.rad, lat=latitude * u.rad, height=0 * u.m
            )
            expected_lst = time.sidereal_time(
                "mean", longitude=longitude * u.rad, model="IAU2006"
            ).to_value(u.rad) % (2 * PI)
            self.assertLess(
                abs(((lst * 2 * PI - expected_lst + PI) % (2 * PI)) - PI) * u.rad,
                0.1 * u.arcsec,
            )

            az_alt = self._probe(
                [["eq2horiz", "h", ra, dec, jd, latitude, longitude]]
            )[0]
            # HADec -> AltAz is Astropy's pure local rotation for mean RA/Dec;
            # starting with FK5 -> AltAz would add apparent/topocentric corrections.
            hour_angle = (expected_lst - ra) * u.rad
            expected = SkyCoord(
                ha=hour_angle,
                dec=dec * u.rad,
                frame=HADec(obstime=time, location=location),
            ).transform_to(AltAz(obstime=time, location=location, pressure=0 * u.hPa))
            actual = SkyCoord(
                az=az_alt[0] * u.rad,
                alt=az_alt[1] * u.rad,
                frame=AltAz(obstime=time, location=location, pressure=0 * u.hPa),
            )
            self.assertLess(_coordinate_separation(actual, expected), 0.1 * u.arcsec)

            inverse = self._probe(
                [["eq2horiz", "e", az_alt[0], az_alt[1], jd, latitude, longitude]]
            )[0]
            original = SkyCoord(ra=ra * u.rad, dec=dec * u.rad, frame=self.fk5)
            reconstructed = SkyCoord(
                ra=_normalized_ra(inverse[0]) * u.rad,
                dec=inverse[1] * u.rad,
                frame=self.fk5,
            )
            self.assertLess(_coordinate_separation(original, reconstructed), 0.01 * u.mas)

        for longitude_deg, latitude_deg, height in [
            (0.0, 0.0, 0.0),
            (45.0, 30.0, 123.0),
            (180.0, -45.0, 2000.0),
            (-122.3, 37.8, 15.0),
        ]:
            longitude, latitude = _radians(longitude_deg), _radians(latitude_deg)
            result = self._probe(
                [["geod2geoc", longitude, latitude, height, "WGS84"]]
            )[0]
            location = EarthLocation.from_geodetic(
                lon=longitude * u.rad,
                lat=latitude * u.rad,
                height=height * u.m,
                ellipsoid="WGS84",
            )
            x, y, z = [value.to_value(u.m) for value in location.to_geocentric()]
            np.testing.assert_allclose(result[3:6], [x, y, z], atol=1e-6)
            expected_spherical = SkyCoord(
                CartesianRepresentation(x=x * u.m, y=y * u.m, z=z * u.m),
                frame=Galactic(),
            ).spherical
            longitude_error = (result[0] - expected_spherical.lon.to_value(u.rad) + PI) % (2 * PI) - PI
            self.assertAlmostEqual(longitude_error, 0.0, delta=1e-12)
            self.assertAlmostEqual(result[1], expected_spherical.lat.to_value(u.rad), delta=1e-12)
            self.assertAlmostEqual(result[2], expected_spherical.distance.to_value(u.m), delta=1e-6)

    def test_angle_unit_and_fold_helpers(self):
        input_value = 1.2
        result = self._probe([["units", input_value]])[0]
        expected = [
            (input_value * u.rad).to_value(u.hourangle),
            (input_value * u.hourangle).to_value(u.rad),
            (input_value * u.rad).to_value(u.deg),
            (input_value * u.deg).to_value(u.rad),
        ]
        np.testing.assert_allclose(result, expected, rtol=0, atol=1e-14)

        fold_cases = [
            (-0.25, 0.0, 2 * PI, 2 * PI - 0.25),
            (2 * PI + 0.25, 0.0, 2 * PI, 0.25),
            (2 * PI, 0.0, 2 * PI, 0.0),
        ]
        actual = self._probe([["fold", x, lower, upper] for x, lower, upper, _ in fold_cases])
        for row, (_, _, _, expected_value) in zip(actual, fold_cases):
            self.assertAlmostEqual(row[0], expected_value, delta=1e-14)


if __name__ == "__main__":
    unittest.main(verbosity=2)
