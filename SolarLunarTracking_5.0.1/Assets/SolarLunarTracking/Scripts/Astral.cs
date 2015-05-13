/*
 * Copyright (c) 2014, Roger Lew (rogerlew.gmail.com)
 * Date: 2/25/2015 - 5/13/2015
 * License: Apache Software License 2.0
 * 
 * The project described was supported by NSF award number IIA-1301792
 * from the NSF Idaho EPSCoR Program and by the National Science Foundation.
 * 
 * 
 * Code ported from:
 * https://pypi.python.org/pypi/astral/0.8.1 under Apache Software License 2.0
 * http://www.mathworks.com/matlabcentral/fileexchange/22992-lunar-azimuth-and-altitude-estimation-algorithm/content/LunarAzEl.m under BSD
 */

using System;
using UnityEngine;
using System.Diagnostics;
using System.Collections;
using System.Collections.Generic;

namespace VTL.SolarLunarTracking
{
    class C // stores Astronomical Constants
    {

        // JDN stands for Julian Day Number
        // Angles here are in degrees

        // 1980 January 0.0 in JDN
        public static double epoch = 2444238.5;

        // Ecliptic longitude of the Sun at epoch 1980.0
        public static double ecliptic_longitude_epoch = 278.833540;

        // Ecliptic longitude of the Sun at perigee
        public static double ecliptic_longitude_perigee = 282.596403;

        // Eccentricity of Earth's orbit
        public static double eccentricity = 0.016718;

        // Semi-major axis of Earth's orbit, in kilometers
        public static double sun_smaxis = 1.49585e8;

        // Sun's angular size, in degrees, at semi-major axis distance
        public static double sun_angular_size_smaxis = 0.533128;

        //// Elements of the Moon's orbit, epoch 1980.0

        // Moon's mean longitude at the epoch
        public static double moon_mean_longitude_epoch = 64.975464;

        // Mean longitude of the perigee at the epoch
        public static double moon_mean_perigee_epoch = 349.383063;

        // Mean longitude of the node at the epoch
        public static double node_mean_longitude_epoch = 151.950429;

        // Inclination of the Moon's orbit
        public static double moon_inclination = 5.145396;

        // Eccentricity of the Moon's orbit
        public static double moon_eccentricity = 0.054900;

        // Moon's angular size at distance a from Earth
        public static double moon_angular_size = 0.5181;

        // Semi-mojor axis of the Moon's orbit, in kilometers
        public static double moon_smaxis = 384401.0;

        // Parallax at a distance a from Earth
        public static double moon_parallax = 0.9507;

        // Synodic month (new Moon to new Moon), in days
        public static double synodic_month = 29.53058868;

        // Base date for E. W. Brown's numbered series of lunations (1923 January 16)
        public static double lunations_base = 2423436.0;

        //// Properties of the Earth
        public static double earth_radius = 6378.16;

    }

    public struct Solar
    {
        public double azimuth;
        public double elevation;
        public double distance;
        public double angular_diameter;
        public DateTime? dawn;
        public DateTime? sunrise;
        public DateTime? sunset;
        public DateTime? dusk;
        public DateTime? solar_noon;

        public double temperature;
        public Color color_temp;
    }

    public struct Lunar
    {
        public double azimuth;
        public double elevation;
        public double distance;
        public double age;
        public double phase;
        public double illuminated;
        public double angular_diameter;
    }

    public class Astral
    {
        #region PROPERTIES
        private TimeSpan utcOffset;
        public TimeSpan utcOffsetProperty
        {
            get { return utcOffset; }
            set { utcOffset = value; }
        }

        private double altitude;
        public double altitudeProperty
        {
            get { return altitude; }
            set { altitude = value; }
        }

        private double longitude;
        public double longitudeProperty
        {
            get { return longitude; }
            set { longitude = value; }
        }

        private double latitude;
        public double latitudeProperty
        {
            get { return latitude; }
            set { latitude = value; }
        }

        private double depression;
        public double depressionProperty
        {
            get { return depression; }
            set { depression = value; }
        }
        #endregion

        static double pi = Math.PI;
        static double D2R = pi / 180;
        static double R2D = 180 / pi;

        static double COLOR_SHIFT_LOW_ELEVATION = -6;
        static double COLOR_SHIFT_HIGH_ELEVATION = 6;
        static double COLOR_SHIFT_DAY_TEMPERATURE = 6500;  // kelvin
        static double COLOR_SHIFT_NIGHT_TEMPERATURE = 3500;

        Dictionary<int, Color> temp2color;

        public Astral(double lng, double lat, double alt,
                      double timeZone, double depression = 6)
        {
            longitudeProperty = lng;
            latitudeProperty = lat;
            altitudeProperty = alt / 1000;

            if (latitudeProperty > 89.8)
                latitudeProperty = 89.8;

            if (latitudeProperty < -89.8)
                latitudeProperty = -89.8;

            utcOffset = new TimeSpan(0, (int)(timeZone * 60), 0);

            depressionProperty = depression;

            loadColorTempData();
        }

        #region STATIC MATH FUNCTIONS
        // Makes porting easier

        public void cart2spher(double x, double y, double z,
                               out double theta, out double phi, out double R)
        {
            double x2 = x * x;
            double y2 = y * y;
            double z2 = z * z;

            theta = atan2(y, x);
            phi = atan2(z, sqrt(x2 + y2));
            R = sqrt(x2 + y2 + z2);
        }

        public void spher2cart(double theta, double phi, double R,
                               out double x, out double y, out double z)
        {
            x = R * cos(phi) * cos(theta);
            y = R * cos(phi) * sin(theta);
            z = R * sin(phi);
        }

        static double degrees(double x)
        {
            return x * R2D;
        }

        static double radians(double x)
        {
            return x * D2R;
        }

        static double fixangle(double a)
        {
            return a - 360.0 * floor(a / 360.0);
        }

        static double sqrt(double x)
        {
            return Math.Sqrt(x);
        }

        static double sum(double[] x)
        {
            double s = 0;
            foreach (double v in x)
            {
                s += v;
            }
            return s;
        }

        static double[,] matrix_transpose(double[,] A)
        {
            int rA = A.GetLength(0);
            int cA = A.GetLength(1);

            double[,] B = new double[cA, rA];

            for (int i = 0; i < cA; i++)
                for (int j = 0; j < rA; j++)
                    B[i, j] = A[j, i];

            return B;
        }

        static double[,] matrix_scalar_multiply(double[,] A, double x)
        {
            int rA = A.GetLength(0);
            int cA = A.GetLength(1);

            double[,] B = new double[rA, cA];

            for (int i = 0; i < rA; i++)
                for (int j = 0; j < cA; j++)
                    B[i, j] = A[i, j] * x;

            return B;
        }

        static double matrix_dot(double[,] A, double[,] B)
        {
            int rA = A.GetLength(0);
            int cA = A.GetLength(1);
            int rB = B.GetLength(0);
            int cB = B.GetLength(1);
            double result = 0;

            System.Diagnostics.Debug.Assert(cA == cB);
            System.Diagnostics.Debug.Assert(rA == rB);

            for (int i = 0; i < rA; i++)
            {
                for (int j = 0; j < cA; j++)
                {
                    result += A[i, j] * B[i, j];
                }
            }
            return result;
        }

        static double[,] matrix_multiply(double[,] A, double[,] B)
        {
            int rA = A.GetLength(0);
            int cA = A.GetLength(1);
            //        int rB = B.GetLength(0);
            int cB = B.GetLength(1);
            double temp = 0;
            double[,] kHasil = new double[rA, cB];

            for (int i = 0; i < rA; i++)
            {
                for (int j = 0; j < cB; j++)
                {
                    temp = 0;
                    for (int k = 0; k < cA; k++)
                    {
                        temp += A[i, k] * B[k, j];
                    }
                    kHasil[i, j] = temp;
                }
            }
            return kHasil;
        }

        static double pow(double x, double y)
        {
            return Math.Pow(x, y);
        }

        static double floor(double x)
        {
            return Math.Floor(x);
        }

        static double ceiling(double x)
        {
            return Math.Ceiling(x);
        }

        static double abs(double x)
        {
            return Math.Abs(x);
        }

        static double cos(double x)
        {
            return Math.Cos(x);
        }

        static double sin(double x)
        {
            return Math.Sin(x);
        }

        static double tan(double x)
        {
            return Math.Tan(x);
        }

        static double atan2(double y, double x)
        {
            return Math.Atan2(y, x);
        }

        static double acos(double x)
        {
            return Math.Acos(x);
        }

        static double asin(double x)
        {
            return Math.Asin(x);
        }

        static double atan(double x)
        {
            return Math.Atan(x);
        }
        #endregion

        #region PUBLIC METHODS
        public DateTime? dawn_utc(DateTime utc)
        {
            /*
            Calculate dawn time in the UTC timezone.
        
            :param date:       Date to calculate for.
            :type date:        datetime.date
            :param latitude:   Latitude - Northern latitudes should be positive
            :type latitude:    float 
            :param longitude:  Longitude - Eastern longitudes should be positive
            :type longitude:   float 
        
            :rtype: date/time in UTC timezone
            */

            try
            {
                return _calc_time(utc, latitudeProperty, longitudeProperty, depressionProperty);
            }
            catch
            {
                // Sun remains below horizon on this day,  at this location.
                return null;
            }
        }

        public DateTime? sunrise_utc(DateTime utc)
        {
            /*
            Calculate sunrise time in the UTC timezone.
        
            :param date:       Date to calculate for.
            :type date:        datetime.date
            :param latitude:   Latitude - Northern latitudes should be positive
            :type latitude:    float 
            :param longitude:  Longitude - Eastern longitudes should be positive
            :type longitude:   float 
        
            :rtype: date/time in UTC timezone
            */

            try
            {
                return _calc_time(utc, latitudeProperty, longitudeProperty, 0.833);
            }
            catch
            {
                // Sun remains below horizon on this day,  at this location.
                return null;
            }

        }

        public DateTime? sunset_utc(DateTime utc)
        {
            /*
            Calculate sunset time in the UTC timezone.
        
            :param date:       Date to calculate for.
            :type date:        datetime.date
            :param latitude:   Latitude - Northern latitudes should be positive
            :type latitude:    float 
            :param longitude:  Longitude - Eastern longitudes should be positive
            :type longitude:   float 
        
            :rtype: date/time in UTC timezone
            */

            try
            {
                return _calc_time(utc, latitudeProperty, longitudeProperty, -0.833);
            }
            catch
            {
                // Sun remains below horizon on this day,  at this location.
                return null;
            }
        }

        public DateTime? dusk_utc(DateTime utc)
        {
            /*
            Calculate dusk time in the UTC timezone.
        
            :param date:       Date to calculate for.
            :type date:        datetime.date
            :param latitude:   Latitude - Northern latitudes should be positive
            :type latitude:    float 
            :param longitude:  Longitude - Eastern longitudes should be positive
            :type longitude:   float 
        
            :rtype: date/time in UTC timezone
            */

            try
            {
                return _calc_time(utc, latitudeProperty, longitudeProperty, depressionProperty);
            }
            catch
            {
                // Sun remains below horizon on this day,  at this location.
                return null;
            }
        }

        public DateTime? solar_noon_utc(DateTime utc)
        {
            /*
            Calculate solar noon time in the UTC timezone.
        
            :param date:       Date to calculate for.
            :type date:        datetime.date
            :param latitude:   Latitude - Northern latitudes should be positive
            :type latitude:    float 
            :param longitude:  Longitude - Eastern longitudes should be positive
            :type longitude:   float 
        
            :rtype: date/time in UTC timezone
            */

            double julianday = _julianday(utc);

            double newt = _jday_to_jcentury(julianday + 0.5 + -longitudeProperty / 360.0);

            double eqtime = _eq_of_time(newt);
            double timeUTC = 720.0 + (-longitude * 4.0) - eqtime;

            timeUTC = timeUTC / 60.0;
            int hour = (int)(timeUTC);
            int minute = (int)((timeUTC - hour) * 60);
            int second = (int)((((timeUTC - hour) * 60) - minute) * 60);

            if (second > 59)
            {
                second -= 60;
                minute += 1;
            }
            else if (second < 0)
            {
                second += 60;
                minute -= 1;
            }

            if (minute > 59)
            {
                minute -= 60;
                hour += 1;
            }
            else if (minute < 0)
            {
                minute += 60;
                hour -= 1;
            }

            if (hour > 23)
            {
                hour -= 24;
                utc += new TimeSpan(24, 0, 0);
            }
            else if (hour < 0)
            {
                hour += 24;
                utc -= new TimeSpan(24, 0, 0);
            }

            DateTime? noon = new DateTime(utc.Year, utc.Month, utc.Day,
                                         hour, minute, second);

            return noon;
        }

        public Solar solar(DateTime date)
        {
            /*
            Calculate the elevation of the sun.
        
            :param dateandtime:       Date/time to calculate for.
            :type dateandtime:        datetime.datetime
            :param latitude:   Latitude - Northern latitudes should be positive
            :type latitude:    float 
            :param longitude:  Longitude - Eastern longitudes should be positive
            :type longitude:   float 
        
            :rtype: Elevation in degrees
            */

            double zone = -utcOffset.Hours;

            DateTime utc_datetime = date - utcOffset;
            double timenow = (double)utc_datetime.Hour +
                             ((double)utc_datetime.Minute / 60.0) +
                             ((double)utc_datetime.Second / 3600.0);

            double JD = _julianday(date);
            double t = _jday_to_jcentury(JD + timenow / 24.0);
            double theta = _sun_declination(t);
            double Etime = _eq_of_time(t);

            double eqtime = Etime;
            double solarDec = theta;   // in degrees

            double solarTimeFix = eqtime - (4.0 * -longitude) + (60 * zone);
            double trueSolarTime = date.Hour * 60.0 + date.Minute +
                date.Second / 60.0 + solarTimeFix;   //    in minutes

            while (trueSolarTime > 1440)
            {
                trueSolarTime = trueSolarTime - 1440;
            }

            double hourangle = trueSolarTime / 4.0 - 180.0;
            //    Thanks to Louis Schwarzmayr for the next line:
            if (hourangle < -180)
                hourangle = hourangle + 360.0;

            double harad = radians(hourangle);

            double csz = sin(radians(latitude)) * sin(radians(solarDec)) +
                         cos(radians(latitude)) * cos(radians(solarDec)) * cos(harad);

            if (csz > 1.0)
                csz = 1.0;

            else if (csz < -1.0)
                csz = -1.0;

            double zenith = degrees(acos(csz));
            double azDenom = (cos(radians(latitude)) * sin(radians(zenith)));

            double solar_az, solar_el, solar_dist, solar_ang_diam;

            if (abs(azDenom) > 0.001)
            {
                double azRad = ((sin(radians(latitude)) * cos(radians(zenith))) -
                                sin(radians(solarDec))) / azDenom;

                if (abs(azRad) > 1.0)
                {
                    if (azRad < 0)
                        azRad = -1.0;
                    else
                        azRad = 1.0;
                }
                solar_az = 180.0 - degrees(acos(azRad));

                if (hourangle > 0.0)
                    solar_az = -solar_az;
            }
            else
            {
                if (latitude > 0.0)
                    solar_az = 180.0;
                else
                    solar_az = 0;
            }
            if (solar_az < 0.0)
                solar_az = solar_az + 360.0;

            double exoatmElevation = 90.0 - zenith;
            double refractionCorrection;

            if (exoatmElevation > 85.0)
            {
                refractionCorrection = 0.0;
            }
            else
            {
                double te = tan(radians(exoatmElevation));
                if (exoatmElevation > 5.0)
                {
                    refractionCorrection = 58.1 / te - 0.07 / (te * te * te) +
                        0.000086 / (te * te * te * te * te);
                }
                else if (exoatmElevation > -0.575)
                {
                    double step1 = (-12.79 + exoatmElevation * 0.711);
                    double step2 = (103.4 + exoatmElevation * (step1));
                    double step3 = (-518.2 + exoatmElevation * (step2));
                    refractionCorrection = 1735.0 + exoatmElevation * (step3);
                }
                else
                {
                    refractionCorrection = -20.774 / te;
                }
                refractionCorrection = refractionCorrection / 3600.0;
            }
            double solarzen = zenith - refractionCorrection;

            solar_el = 90.0 - solarzen;

            double F = solar_orbital_distance_factor(date);

            // Distance to Sun in km
            solar_dist = C.sun_smaxis / F;
            solar_ang_diam = F * C.sun_angular_size_smaxis;

            Solar _solar = new Solar();
            _solar.azimuth = solar_az;
            _solar.elevation = solar_el;
            _solar.distance = solar_dist;
            _solar.angular_diameter = solar_ang_diam;

            _solar.dawn = dawn_utc(utc_datetime) + utcOffset;
            _solar.sunrise = sunrise_utc(utc_datetime) + utcOffset;
            _solar.sunset = sunset_utc(utc_datetime) + utcOffset;
            _solar.dusk = dusk_utc(utc_datetime) + utcOffset;
            _solar.solar_noon = solar_noon_utc(utc_datetime) + utcOffset;

            _solar.temperature = solar_temperature(solar_el);
            _solar.color_temp = temp2Color(_solar.temperature);

            return _solar;
        }

        DateTime utc(DateTime local)
        {
            return local - utcOffset;
        }

        DateTime local(DateTime utc)
        {
            return utc + utcOffset;
        }

        public Lunar lunar(DateTime date)
        {
            DateTime utc_datetime = date - utcOffset;

            // Convert Universal Time to Ephemeris Time
            double jd = _juliandate(utc_datetime);

            // Find the Day Number
            double d = jd - 2451543.5;


            // Keplerian Elements of the Moon
            // This will also account for the Sun's perturbation
            double N = 125.1228 - 0.0529538083 * d;             // (Long asc. node deg)
            double i = 5.1454;                                  // (Inclination deg)
            double w = 318.0634 + 0.1643573223 * d;               // (Arg. of perigee deg)
            double a = 60.2666;                                // (Mean distance (Earth's Equitorial Radii)
            double e = 0.054900;                                // (Eccentricity)
            double M = (115.3654 + 13.0649929509 * d) % 360;    // (Mean anomaly deg)

            double LMoon = (N + w + M) % 360;                   // (Moon's mean longitude deg)
            double FMoon = (LMoon - N) % 360;                   // (Moon's argument of latitude)



            // Keplerian Elements of the Sun
            double wSun = (282.9404 + 4.70935E-5 * d) % 360;    // (longitude of perihelion)
            double MSun = (356.0470 + 0.9856002585 * d) % 360;  // (Sun mean anomaly)
            double LSun = (wSun + MSun) % 360;                  // (Sun's mean longitude)

            double DMoon = LMoon - LSun;                       // (Moon's mean elongation)  


            // Calculate Lunar perturbations in Longitude
            double[] LunarPLon = new double[] 
                { -1.274*sin(radians(M - 2.0 * DMoon)),
                   0.658*sin(radians(2.0 * DMoon)),
                  -0.186*sin(radians(MSun)),
                  -0.059*sin(radians(2.0 * M - 2.0 * DMoon)),
                  -0.057*sin(radians(M - 2.0 * DMoon + MSun)),
                   0.053*sin(radians(M + 2.0 * DMoon)),
                   0.046*sin(radians(2.0 * DMoon - MSun)),
                   0.041*sin(radians(M - MSun)),
                  -0.035*sin(radians(DMoon)),           
                  -0.031*sin(radians(M + MSun)),
                  -0.015*sin(radians(2.0 * FMoon-2.0 * DMoon)),
                   0.011*sin(radians(M - 4.0 * DMoon)) };

            // Calculate Lunar perturbations in Latitude 
            double[] LunarPLat = new double[]
                { -0.173 * sin(radians(FMoon - 2.0 * DMoon)),
                  -0.055 * sin(radians(M - FMoon - 2 * DMoon)),
                  -0.046 * sin(radians(M + FMoon - 2.0 * DMoon)),
                  +0.033 * sin(radians(FMoon + 2.0 * DMoon)),
                  +0.017 * sin(radians(2.0 * M + FMoon)) };

            // Calculate perturbations in Distance
            double[] LunarPDist = new double[]
                { -0.58*cos(radians(M - 2.0 * DMoon)),
                  -0.46*cos(radians(2*DMoon)) };

            // Compute E, the eccentric anomaly

            // E0 is the eccentric anomaly approximation estimate 
            // (this will initially have a relativly high error)
            double radM = radians(M);
            double E0 = M + (180.0 / pi) * e * sin(radM) * (1 + e * cos(radM));

            // Compute E1 and set it to E0 until the E1 == E0
            double radE0 = radians(E0);
            double E1 = E0 - (E0 - (180 / pi) * e * sin(radE0) - M) / (1 - e * cos(radE0));
            while (E1 - E0 > .000005)
            {
                E0 = E1;
                radE0 = radians(E0);
                E1 = E0 - (E0 - (180 / pi) * e * sin(radE0) - M) / (1 - e * cos(radE0));
            }

            double E = E1;

            // Compute rectangular coordinates (x,y) in the plane of the lunar orbit
            double radE = radians(E);
            double x = a * (cos(radE) - e);
            double y = a * sqrt(1 - e * e) * sin(radE);

            // Convert this to distance and true anomaly
            double r = sqrt(x * x + y * y);
            double v = degrees(atan2(radians(y), radians(x)));

            // Compute moon's position in ecliptic coordinates
            double radN = radians(N);
            double rad_vw = radians(v + w);
            double rad_i = radians(i);

            double xeclip = r * (cos(radN) * cos(rad_vw) - sin(radN) * sin(rad_vw) * cos(rad_i));
            double yeclip = r * (sin(radN) * cos(rad_vw) + cos(radN) * sin((rad_vw)) * cos(rad_i));
            double zeclip = r * sin(rad_vw) * sin(rad_i);

            // Add the calculated lunar perturbation terms to increase model fidelity
            double eLon, eLat, eDist;
            cart2spher(xeclip, yeclip, zeclip, out eLon, out eLat, out eDist);

            spher2cart(eLon + radians(sum(LunarPLon)),
                       eLat + radians(sum(LunarPLat)),
                       eDist + sum(LunarPDist),
                       out xeclip, out yeclip, out zeclip);

            // Convert the latitude and longitude to right ascension RA and declination
            // delta
            double T = (jd - 2451545.0) / 36525.0;

            // Generate a rotation matrix for ecliptic to equitorial
            // RotM=rotm_coo('E',jd)
            // See rotm_coo.m for obl and rotational matrix transformation
            double T2 = T * T;
            double Obl = 23.439291 - 0.0130042 * T - 0.00000016 * T2 + 0.000000504 * T2 * T;
            Obl = radians(Obl);

            double[,] RotM = new double[,] 
                  { {1, 0, 0}, 
                    {0, cos(Obl), sin(Obl)}, 
                    {0, -sin(Obl), cos(Obl)} };

            RotM = matrix_transpose(RotM);

            // Apply the rotational matrix to the ecliptic rectangular coordinates
            // Also, convert units to km instead of earth equatorial radii
            double[,] eclipticCoords = new double[,] { { xeclip }, { yeclip }, { zeclip } };
            double[,] sol2 = matrix_multiply(RotM, eclipticCoords);
            sol2 = matrix_scalar_multiply(sol2, C.earth_radius);
            double[] sol = new double[] { sol2[0, 0], sol2[1, 0], sol2[2, 0] };

            // Find the equatorial rectangular coordinates of the location specified
            double rad_lng = radians(longitude);
            double rad_lat = radians(latitude);
            double xel, yel, zel;
            spher2cart(rad_lng, rad_lat, altitude + C.earth_radius,
                       out xel, out yel, out zel);

            // Find the equatorial rectangular coordinates of the location @ sea level
            double xsl, ysl, zsl;
            spher2cart(rad_lng, rad_lat, C.earth_radius,
                       out xsl, out ysl, out zsl);

            // Find the Angle Between sea level coordinate vector and the moon vector
            double[,] A1 = new double[,] { { xsl, ysl, zsl } };
            double[,] B1 = new double[,] { { sol[0] - xsl, sol[1] - ysl, sol[2] - zsl } };
            double theta1 = matrix_dot(A1, B1);
            theta1 /= sqrt(pow(xsl, 2) + pow(ysl, 2) + pow(zsl, 2)) *
                      sqrt(pow(sol[0] - xsl, 2) + pow(sol[1] - ysl, 2) + pow(sol[2] - zsl, 2));
            theta1 = 180 - degrees(acos(theta1));

            // Find the Angle Between the same coordinates but at the specified elevation
            double[,] A2 = new double[,] { { xel, yel, zel } };
            double[,] B2 = new double[,] { { sol[0] - xel, sol[1] - yel, sol[2] - zel } };
            double theta2 = matrix_dot(A2, B2);
            theta2 /= sqrt(pow(xel, 2) + pow(yel, 2) + pow(zel, 2)) *
                      sqrt(pow(sol[0] - xel, 2) + pow(sol[1] - yel, 2) + pow(sol[2] - zel, 2));
            theta2 = 180 - degrees(acos(theta2));

            // Find the Difference Between the two angles (+|-) is important
            double thetaDiff = theta2 - theta1;

            // Equatorial to horizon coordinate transformation
            double RA, delta, _phi;
            cart2spher(sol[0], sol[1], sol[2], out RA, out delta, out _phi);
            delta = degrees(delta);
            RA = degrees(RA);

            // Following the RA DEC to Az Alt conversion sequence explained here:
            // http://www.stargazing.net/kepler/altaz.html

            // Find the J2000 value
            double J2000 = jd - 2451545.0;
            double UTH = utc_datetime.Hour + utc_datetime.Minute / 60.0 + utc_datetime.Second / 3600.0;

            //Calculate local siderial time
            double LST = (100.46 + 0.985647 * J2000 + longitude + 15 * UTH) % 360;

            // Replace RA with hour angle HA
            double HA = LST - RA;

            // Find the h and AZ at the current LST
            double rad_delta = radians(delta);
            double rad_HA = radians(HA);
            double rad_h = asin(sin(rad_delta) * sin(rad_lat) +
                           cos(rad_delta) * cos(rad_lat) * cos(rad_HA));
            double lunar_el = degrees(rad_h);
            double rad_Az = acos((sin(rad_delta) -
                            (sin(rad_h)) * sin(rad_lat)) / (cos(rad_h) * cos(rad_lat)));
            double lunar_az = degrees(rad_Az);

            // Add in the angle offset due to the specified site elevation
            lunar_el += thetaDiff;

            if (sin(rad_HA) >= 0)
                lunar_az = 360 - lunar_az;

            // Apply Paralax Correction if we are still on earth
            if (altitude < 100)
            {
                double horParal = 8.794 / (r * 6379.14 / 149.59787e6);
                double p = degrees(asin(cos(radians(lunar_el)) * sin(radians(horParal / 3600))));
                lunar_el -= p;
            }

            Lunar _lunar = new Lunar();
            _lunar.azimuth = lunar_az;
            _lunar.elevation = lunar_el;

            lunar_other(date, ref _lunar);
            return _lunar;
        }

        static double kepler(double m, double ecc, double epsilon = 1e-6)
        {
            m = radians(m);
            double e = m;
            double delta;
            while (true)
            {
                delta = e - ecc * sin(e) - m;
                e -= delta / (1.0 - ecc * cos(e));

                if (abs(delta) <= epsilon)
                    break;
            }
            return e;
        }

        static double solar_orbital_distance_factor(DateTime date)
        {
            // Used to calculate solar distance and viewing angle
            // uses a different epoch from solar. Their is some redundancy in calculation
            // but it works, so I think it is better to just leave it alone.

            double day = _juliandate(date) - C.epoch;

            // Mean anomaly of the Sun
            double N = fixangle((360 / 365.2422) * day);

            // Convert from perigee coordinates to epoch 1980
            double M = fixangle(N + C.ecliptic_longitude_epoch - C.ecliptic_longitude_perigee);

            // Solve Kepler's equation
            double Ec = kepler(M, C.eccentricity);
            Ec = sqrt((1 + C.eccentricity) / (1 - C.eccentricity)) * tan(Ec / 2.0);

            // True anomaly
            Ec = 2 * degrees(atan(Ec));

            //        // Suns's geometric ecliptic longuitude
            //        double lambda_sun = fixangle(Ec + C.ecliptic_longitude_perigee);

            // Orbital distance factor
            double F = ((1 + C.eccentricity * cos(radians(Ec))) / (1 - pow(C.eccentricity, 2)));

            return F;
        }

        void lunar_other(DateTime date, ref Lunar _lunar)
        {
            // Calculates lunar distance, phase, and angular diameter
            // uses a different epoch from lunar. Their is some redundancy in calculation
            // but it works, so I think it is better to just leave it alone.


            double day = _juliandate(date) - C.epoch;

            // Mean anomaly of the Sun
            double N = fixangle((360 / 365.2422) * day);

            // Convert from perigee coordinates to epoch 1980
            double M = fixangle(N + C.ecliptic_longitude_epoch - C.ecliptic_longitude_perigee);

            // Solve Kepler's equation
            double Ec = kepler(M, C.eccentricity);
            Ec = sqrt((1 + C.eccentricity) / (1 - C.eccentricity)) * tan(Ec / 2.0);

            // True anomaly
            Ec = 2 * degrees(atan(Ec));

            // Suns's geometric ecliptic longuitude
            double lambda_sun = fixangle(Ec + C.ecliptic_longitude_perigee);

            // Orbital distance factor
            //        double F = ((1 + C.eccentricity * cos(radians(Ec))) / (1 - pow(C.eccentricity, 2)));

            //       // Distance to Sun in km
            //       double sun_dist = C.sun_smaxis / F;
            //       double sun_angular_diameter = F * C.sun_angular_size_smaxis;

            ////////////////
            //
            // Calculation of the Moon's position

            // Moon's mean longitude
            double moon_longitude = fixangle(13.1763966 * day + C.moon_mean_longitude_epoch);

            // Moon's mean anomaly
            double MM = fixangle(moon_longitude - 0.1114041 * day - C.moon_mean_perigee_epoch);

            // Moon's ascending node mean longitude
            // MN = fixangle(c.node_mean_longitude_epoch - 0.0529539 * day)

            double evection = 1.2739 * sin(radians(2 * (moon_longitude - lambda_sun) - MM));

            // Annual equation
            double annual_eq = 0.1858 * sin(radians(M));

            // Correction term
            double A3 = 0.37 * sin(radians(M));
            double MmP = MM + evection - annual_eq - A3;

            // Correction for the equation of the centre
            double mEc = 6.2886 * sin(radians(MmP));

            // Another correction term
            double A4 = 0.214 * sin(radians(2 * MmP));

            // Corrected longitude
            double lP = moon_longitude + evection + mEc - annual_eq + A4;

            // Variation
            double variation = 0.6583 * sin(radians(2 * (lP - lambda_sun)));

            // True longitude
            double lPP = lP + variation;

            //
            // Calculation of the Moon's inclination
            // unused for phase calculation.

            // Corrected longitude of the node
            // NP = MN - 0.16 * sin(radians(M))

            // Y inclination coordinate
            // y = sin(radians(lPP - NP)) * cos(radians(c.moon_inclination))

            // X inclination coordinate
            // x = cos(radians(lPP - NP))

            // Ecliptic longitude (unused?)
            // lambda_moon = todeg(atan2(y,x)) + NP

            // Ecliptic latitude (unused?)
            // BetaM = todeg(asin(sin(radians(lPP - NP)) * sin(radians(c.moon_inclination))))

            //////////////
            //
            // Calculation of the phase of the Moon

            // Age of the Moon, in degrees
            double moon_age = lPP - lambda_sun;

            // Phase of the Moon
            double moon_phase = (1 - cos(radians(moon_age))) / 2.0;

            // Calculate distance of Moon from the centre of the Earth
            double moon_dist = (C.moon_smaxis * (1 - pow(C.moon_eccentricity, 2))) /
                               (1 + C.moon_eccentricity * cos(radians(MmP + mEc)));

            // Calculate Moon's angular diameter
            double moon_diam_frac = moon_dist / C.moon_smaxis;
            double moon_angular_diameter = C.moon_angular_size / moon_diam_frac;

            _lunar.age = C.synodic_month * fixangle(moon_age) / 360.0;
            _lunar.phase = fixangle(moon_age) / 360.0;
            _lunar.illuminated = moon_phase;
            _lunar.distance = moon_dist;
            _lunar.angular_diameter = moon_angular_diameter;

        }
        #endregion

        #region STATIC METHODS
        static double _proper_angle(double value)
        {
            if (value > 0.0)
            {
                value /= 360.0;
                return (value - floor(value)) * 360;
            }
            else
            {
                double tmp = ceiling(abs(value / 360.0));
                return value + tmp * 360.0;
            }
        }

        static double _julianday(DateTime date)
        {
            // Convert Universal Time to Ephemeris Time

            double day = date.Day;
            double month = date.Month;
            double year = date.Year;

            //if timezone is not None:
            //    offset = timezone.localize(datetime.datetime(year, month, day)).utcoffset()
            //    offset = offset.total_seconds() / 1440.0
            //    day += offset + 0.5

            if (month <= 2)
            {
                year = year - 1;
                month = month + 12;
            }

            double A = floor(year / 100.0);
            double B = 2.0 - A + floor(A / 4.0);

            double jd = floor(365.25 * (year + 4716)) +
                        floor(30.6001 * (month + 1)) +
                        day - 1524.5;

            if (jd > 2299160.4999999)
                jd += B;

            return jd;
        }

        static double _juliandate(DateTime utc)
        {

            double y = utc.Year;
            double m = utc.Month;
            double d = utc.Day;
            double hr = utc.Hour;
            double min = utc.Minute;
            double sec = utc.Second;

            // Address differences between python/C# and C time conventions
            //       C:                Python datetime
            // 0 <= mon  <= 11        1 <= month <= 12
            //

            // C code to get the julian date of the start of the day */
            // takes as input 1900+ptm->tm_year, ptm->tm_mon+1, ptm->tm_mday
            // So we can use just (year, month, mday)

            int mterm = (int)((m - 14) / 12);
            int aterm = (int)((1461 * (y + 4800 + mterm)) / 4);
            int bterm = (int)((367 * (m - 2 - 12 * mterm)) / 12);
            int cterm = (int)((3 * (int)((y + 4900 + mterm) / 100)) / 4);

            double j = aterm + bterm - cterm + d;
            j -= 32075;

            // offset to start of day
            j -= 0.5;

            // Apply the time
            double jd = j + (hr + (min + (sec / 60.0)) / 60.0) / 24.0;

            return jd;

        }

        static double _jday_to_jcentury(double julianday)
        {
            return (julianday - 2451545.0) / 36525.0;
        }

        static double _jcentury_to_jday(double juliancentury)
        {
            return (juliancentury * 36525.0) + 2451545.0;
        }

        static double _mean_obliquity_of_ecliptic(double juliancentury)
        {
            double seconds = 21.448 - juliancentury *
                (46.815 + juliancentury * (0.00059 - juliancentury * (0.001813)));
            return 23.0 + (26.0 + (seconds / 60.0)) / 60.0;
        }

        static double _obliquity_correction(double juliancentury)
        {
            double e0 = _mean_obliquity_of_ecliptic(juliancentury);

            double omega = 125.04 - 1934.136 * juliancentury;
            return e0 + 0.00256 * cos(radians(omega));
        }

        static double _geom_mean_long_sun(double juliancentury)
        {
            double l0 = 280.46646 +
                juliancentury * (36000.76983 + 0.0003032 * juliancentury);
            return l0 % 360.0;
        }

        static double _eccentrilocation_earth_orbit(double juliancentury)
        {
            return 0.016708634 -
                juliancentury * (0.000042037 + 0.0000001267 * juliancentury);
        }

        static double _geom_mean_anomaly_sun(double juliancentury)
        {
            return 357.52911 +
                juliancentury * (35999.05029 - 0.0001537 * juliancentury);
        }

        static double _eq_of_time(double juliancentury)
        {
            double epsilon = _obliquity_correction(juliancentury);
            double l0 = _geom_mean_long_sun(juliancentury);
            double e = _eccentrilocation_earth_orbit(juliancentury);
            double m = _geom_mean_anomaly_sun(juliancentury);

            double y = tan(radians(epsilon) / 2.0);
            y = y * y;

            double sin2l0 = sin(2.0 * radians(l0));
            double sinm = sin(radians(m));
            double cos2l0 = cos(2.0 * radians(l0));
            double sin4l0 = sin(4.0 * radians(l0));
            double sin2m = sin(2.0 * radians(m));

            double Etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0 -
                    0.5 * y * y * sin4l0 - 1.25 * e * e * sin2m;

            return degrees(Etime) * 4.0;
        }

        static double _sun_eq_of_center(double juliancentury)
        {
            double m = _geom_mean_anomaly_sun(juliancentury);

            double mrad = radians(m);
            double sinm = sin(mrad);
            double sin2m = sin(mrad + mrad);
            double sin3m = sin(mrad + mrad + mrad);

            double c = sinm * (1.914602 - juliancentury *
                       (0.004817 + 0.000014 * juliancentury)) +
                       sin2m * (0.019993 - 0.000101 * juliancentury) + sin3m * 0.000289;

            return c;
        }

        static double _sun_true_long(double juliancentury)
        {
            double l0 = _geom_mean_long_sun(juliancentury);
            double c = _sun_eq_of_center(juliancentury);

            return l0 + c;
        }

        static double _sun_apparent_long(double juliancentury)
        {
            double O = _sun_true_long(juliancentury);

            double omega = 125.04 - 1934.136 * juliancentury;
            return O - 0.00569 - 0.00478 * sin(radians(omega));
        }

        static double _sun_declination(double juliancentury)
        {
            double e = _obliquity_correction(juliancentury);
            double lambd = _sun_apparent_long(juliancentury);

            double sint = sin(radians(e)) * sin(radians(lambd));
            return degrees(asin(sint));
        }

        static double _sun_rad_vector(double juliancentury)
        {
            double v = _sun_true_anomoly(juliancentury);
            double e = _eccentrilocation_earth_orbit(juliancentury);

            return (1.000001018 * (1 - e * e)) / (1 + e * cos(radians(v)));
        }

        static double _sun_rt_ascension(double juliancentury)
        {
            double e = _obliquity_correction(juliancentury);
            double lambd = _sun_apparent_long(juliancentury);

            double tananum = (cos(radians(e)) * sin(radians(lambd)));
            double tanadenom = (cos(radians(lambd)));

            return degrees(atan2(tananum, tanadenom));
        }

        static double _sun_true_anomoly(double juliancentury)
        {
            double m = _geom_mean_anomaly_sun(juliancentury);
            double c = _sun_eq_of_center(juliancentury);

            return m + c;
        }

        static double _hour_angle(double latitude, double solar_dec, double solar_depression)
        {
            double latRad = radians(latitude);
            double sdRad = radians(solar_dec);

            double HA = (acos(cos(radians(90 + solar_depression)) /
                        (cos(latRad) * cos(sdRad)) - tan(latRad) * tan(sdRad)));

            return HA;
        }

        public static DateTime _calc_time(DateTime utc, double latitude, double longitude, double depression)
        {
            double julianday = _julianday(utc);

            double t = _jday_to_jcentury(julianday);
            double eqtime = _eq_of_time(t);
            double solarDec = _sun_declination(t);

            double hourangle = -_hour_angle(latitude, solarDec, 0.833);

            double delta = -longitude - degrees(hourangle);
            double timeDiff = 4.0 * delta;
            double timeUTC = 720.0 + timeDiff - eqtime;

            double newt = _jday_to_jcentury(_jcentury_to_jday(t) +
                timeUTC / 1440.0);
            eqtime = _eq_of_time(newt);
            solarDec = _sun_declination(newt);

            if (depression < 0)
            {
                depression = abs(depression);
                hourangle = -_hour_angle(latitude, solarDec, depression);
            }
            else
            {
                hourangle = _hour_angle(latitude, solarDec, depression);
            }

            delta = -longitude - degrees(hourangle);
            timeDiff = 4 * delta;
            timeUTC = 720 + timeDiff - eqtime;

            timeUTC = timeUTC / 60.0;
            int hour = (int)timeUTC;
            int minute = (int)((timeUTC - hour) * 60);
            int second = (int)((((timeUTC - hour) * 60) - minute) * 60);

            if (second > 59)
            {
                second -= 60;
                minute += 1;
            }
            else if (second < 0)
            {
                second += 60;
                minute -= 1;
            }

            if (minute > 59)
            {
                minute -= 60;
                hour += 1;
            }
            else if (minute < 0)
            {
                minute += 60;
                hour -= 1;
            }

            if (hour > 23)
            {
                hour -= 24;
                utc += new TimeSpan(1, 0, 0, 0);
            }
            else if (hour < 0)
                hour += 24;
            utc -= new TimeSpan(1, 0, 0, 0);

            DateTime dt = new DateTime(utc.Year, utc.Month, utc.Day,
                                       hour, minute, second);
            return dt;
        }
        #endregion

        public double solar_temperature(double elevation)
        {
            double temp;
            if (elevation < COLOR_SHIFT_LOW_ELEVATION)
            {
                temp = COLOR_SHIFT_NIGHT_TEMPERATURE;
            }
            else if (elevation < COLOR_SHIFT_HIGH_ELEVATION)
            {
                /* Transition period: interpolate */
                double a = (elevation - COLOR_SHIFT_LOW_ELEVATION) /
                           (COLOR_SHIFT_HIGH_ELEVATION - COLOR_SHIFT_LOW_ELEVATION);

                temp = a * (COLOR_SHIFT_DAY_TEMPERATURE - COLOR_SHIFT_NIGHT_TEMPERATURE) + COLOR_SHIFT_NIGHT_TEMPERATURE;
            }
            else
            {
                temp = COLOR_SHIFT_DAY_TEMPERATURE;
            }
            return temp;
        }

        Color temp2Color(double temp)
        {
            int intTemp = (int)Math.Round(temp / 100) * 100;
            return temp2color[intTemp];
        }

        void loadColorTempData()
        {
            temp2color = new Dictionary<int, Color>();
            temp2color.Add(1000, new Color(1.00000000f, 0.18172716f, 0.00000000f));
            temp2color.Add(1100, new Color(1.00000000f, 0.25503671f, 0.00000000f));
            temp2color.Add(1200, new Color(1.00000000f, 0.30942099f, 0.00000000f));
            temp2color.Add(1300, new Color(1.00000000f, 0.35357379f, 0.00000000f));
            temp2color.Add(1400, new Color(1.00000000f, 0.39091524f, 0.00000000f));
            temp2color.Add(1500, new Color(1.00000000f, 0.42322816f, 0.00000000f));
            temp2color.Add(1600, new Color(1.00000000f, 0.45159884f, 0.00000000f));
            temp2color.Add(1700, new Color(1.00000000f, 0.47675916f, 0.00000000f));
            temp2color.Add(1800, new Color(1.00000000f, 0.49923747f, 0.00000000f));
            temp2color.Add(1900, new Color(1.00000000f, 0.51943421f, 0.00000000f));
            temp2color.Add(2000, new Color(1.00000000f, 0.54360078f, 0.08679949f));
            temp2color.Add(2100, new Color(1.00000000f, 0.56618736f, 0.14065513f));
            temp2color.Add(2200, new Color(1.00000000f, 0.58734976f, 0.18362641f));
            temp2color.Add(2300, new Color(1.00000000f, 0.60724493f, 0.22137978f));
            temp2color.Add(2400, new Color(1.00000000f, 0.62600248f, 0.25591950f));
            temp2color.Add(2500, new Color(1.00000000f, 0.64373109f, 0.28819679f));
            temp2color.Add(2600, new Color(1.00000000f, 0.66052319f, 0.31873863f));
            temp2color.Add(2700, new Color(1.00000000f, 0.67645822f, 0.34786758f));
            temp2color.Add(2800, new Color(1.00000000f, 0.69160518f, 0.37579588f));
            temp2color.Add(2900, new Color(1.00000000f, 0.70602449f, 0.40267128f));
            temp2color.Add(3000, new Color(1.00000000f, 0.71976951f, 0.42860152f));
            temp2color.Add(3100, new Color(1.00000000f, 0.73288760f, 0.45366838f));
            temp2color.Add(3200, new Color(1.00000000f, 0.74542112f, 0.47793608f));
            temp2color.Add(3300, new Color(1.00000000f, 0.75740814f, 0.50145662f));
            temp2color.Add(3400, new Color(1.00000000f, 0.76888303f, 0.52427322f));
            temp2color.Add(3500, new Color(1.00000000f, 0.77987699f, 0.54642268f));
            temp2color.Add(3600, new Color(1.00000000f, 0.79041843f, 0.56793692f));
            temp2color.Add(3700, new Color(1.00000000f, 0.80053332f, 0.58884417f));
            temp2color.Add(3800, new Color(1.00000000f, 0.81024551f, 0.60916971f));
            temp2color.Add(3900, new Color(1.00000000f, 0.81957693f, 0.62893653f));
            temp2color.Add(4000, new Color(1.00000000f, 0.82854786f, 0.64816570f));
            temp2color.Add(4100, new Color(1.00000000f, 0.83717703f, 0.66687674f));
            temp2color.Add(4200, new Color(1.00000000f, 0.84548188f, 0.68508786f));
            temp2color.Add(4300, new Color(1.00000000f, 0.85347859f, 0.70281616f));
            temp2color.Add(4400, new Color(1.00000000f, 0.86118227f, 0.72007777f));
            temp2color.Add(4500, new Color(1.00000000f, 0.86860704f, 0.73688797f));
            temp2color.Add(4600, new Color(1.00000000f, 0.87576611f, 0.75326132f));
            temp2color.Add(4700, new Color(1.00000000f, 0.88267187f, 0.76921169f));
            temp2color.Add(4800, new Color(1.00000000f, 0.88933596f, 0.78475236f));
            temp2color.Add(4900, new Color(1.00000000f, 0.89576933f, 0.79989606f));
            temp2color.Add(5000, new Color(1.00000000f, 0.90198230f, 0.81465502f));
            temp2color.Add(5100, new Color(1.00000000f, 0.90963069f, 0.82838210f));
            temp2color.Add(5200, new Color(1.00000000f, 0.91710889f, 0.84190889f));
            temp2color.Add(5300, new Color(1.00000000f, 0.92441842f, 0.85523742f));
            temp2color.Add(5400, new Color(1.00000000f, 0.93156127f, 0.86836903f));
            temp2color.Add(5500, new Color(1.00000000f, 0.93853986f, 0.88130458f));
            temp2color.Add(5600, new Color(1.00000000f, 0.94535695f, 0.89404470f));
            temp2color.Add(5700, new Color(1.00000000f, 0.95201559f, 0.90658983f));
            temp2color.Add(5800, new Color(1.00000000f, 0.95851906f, 0.91894041f));
            temp2color.Add(5900, new Color(1.00000000f, 0.96487079f, 0.93109690f));
            temp2color.Add(6000, new Color(1.00000000f, 0.97107439f, 0.94305985f));
            temp2color.Add(6100, new Color(1.00000000f, 0.97713351f, 0.95482993f));
            temp2color.Add(6200, new Color(1.00000000f, 0.98305189f, 0.96640795f));
            temp2color.Add(6300, new Color(1.00000000f, 0.98883326f, 0.97779486f));
            temp2color.Add(6400, new Color(1.00000000f, 0.99448139f, 0.98899179f));
            temp2color.Add(6500, new Color(1.00000000f, 1.00000000f, 1.00000000f));
            temp2color.Add(6600, new Color(0.98947904f, 0.99348723f, 1.00000000f));
            temp2color.Add(6700, new Color(0.97940448f, 0.98722715f, 1.00000000f));
            temp2color.Add(6800, new Color(0.96975025f, 0.98120637f, 1.00000000f));
            temp2color.Add(6900, new Color(0.96049223f, 0.97541240f, 1.00000000f));
            temp2color.Add(7000, new Color(0.95160805f, 0.96983355f, 1.00000000f));
            temp2color.Add(7100, new Color(0.94303638f, 0.96443333f, 1.00000000f));
            temp2color.Add(7200, new Color(0.93480451f, 0.95923080f, 1.00000000f));
            temp2color.Add(7300, new Color(0.92689056f, 0.95421394f, 1.00000000f));
            temp2color.Add(7400, new Color(0.91927697f, 0.94937330f, 1.00000000f));
            temp2color.Add(7500, new Color(0.91194747f, 0.94470005f, 1.00000000f));
            temp2color.Add(7600, new Color(0.90488690f, 0.94018594f, 1.00000000f));
            temp2color.Add(7700, new Color(0.89808115f, 0.93582323f, 1.00000000f));
            temp2color.Add(7800, new Color(0.89151710f, 0.93160469f, 1.00000000f));
            temp2color.Add(7900, new Color(0.88518247f, 0.92752354f, 1.00000000f));
            temp2color.Add(8000, new Color(0.87906581f, 0.92357340f, 1.00000000f));
            temp2color.Add(8100, new Color(0.87315640f, 0.91974827f, 1.00000000f));
            temp2color.Add(8200, new Color(0.86744421f, 0.91604254f, 1.00000000f));
            temp2color.Add(8300, new Color(0.86191983f, 0.91245088f, 1.00000000f));
            temp2color.Add(8400, new Color(0.85657444f, 0.90896831f, 1.00000000f));
            temp2color.Add(8500, new Color(0.85139976f, 0.90559011f, 1.00000000f));
            temp2color.Add(8600, new Color(0.84638799f, 0.90231183f, 1.00000000f));
            temp2color.Add(8700, new Color(0.84153180f, 0.89912926f, 1.00000000f));
            temp2color.Add(8800, new Color(0.83682430f, 0.89603843f, 1.00000000f));
            temp2color.Add(8900, new Color(0.83225897f, 0.89303558f, 1.00000000f));
            temp2color.Add(9000, new Color(0.82782969f, 0.89011714f, 1.00000000f));
            temp2color.Add(9100, new Color(0.82353066f, 0.88727974f, 1.00000000f));
            temp2color.Add(9200, new Color(0.81935641f, 0.88452017f, 1.00000000f));
            temp2color.Add(9300, new Color(0.81530175f, 0.88183541f, 1.00000000f));
            temp2color.Add(9400, new Color(0.81136180f, 0.87922257f, 1.00000000f));
            temp2color.Add(9500, new Color(0.80753191f, 0.87667891f, 1.00000000f));
            temp2color.Add(9600, new Color(0.80380769f, 0.87420182f, 1.00000000f));
            temp2color.Add(9700, new Color(0.80018497f, 0.87178882f, 1.00000000f));
            temp2color.Add(9800, new Color(0.79665980f, 0.86943756f, 1.00000000f));
            temp2color.Add(9900, new Color(0.79322843f, 0.86714579f, 1.00000000f));
            temp2color.Add(10000, new Color(0.78988728f, 0.86491137f, 1.00000000f));
            temp2color.Add(10100, new Color(0.78663296f, 0.86273225f, 1.00000000f));
            temp2color.Add(10200, new Color(0.78346225f, 0.86060650f, 1.00000000f));
            temp2color.Add(10300, new Color(0.78037207f, 0.85853224f, 1.00000000f));
            temp2color.Add(10400, new Color(0.77735950f, 0.85650771f, 1.00000000f));
            temp2color.Add(10500, new Color(0.77442176f, 0.85453121f, 1.00000000f));
            temp2color.Add(10600, new Color(0.77155617f, 0.85260112f, 1.00000000f));
            temp2color.Add(10700, new Color(0.76876022f, 0.85071588f, 1.00000000f));
            temp2color.Add(10800, new Color(0.76603147f, 0.84887402f, 1.00000000f));
            temp2color.Add(10900, new Color(0.76336762f, 0.84707411f, 1.00000000f));
            temp2color.Add(11000, new Color(0.76076645f, 0.84531479f, 1.00000000f));
            temp2color.Add(11100, new Color(0.75822586f, 0.84359476f, 1.00000000f));
            temp2color.Add(11200, new Color(0.75574383f, 0.84191277f, 1.00000000f));
            temp2color.Add(11300, new Color(0.75331843f, 0.84026762f, 1.00000000f));
            temp2color.Add(11400, new Color(0.75094780f, 0.83865816f, 1.00000000f));
            temp2color.Add(11500, new Color(0.74863017f, 0.83708329f, 1.00000000f));
            temp2color.Add(11600, new Color(0.74636386f, 0.83554194f, 1.00000000f));
            temp2color.Add(11700, new Color(0.74414722f, 0.83403311f, 1.00000000f));
            temp2color.Add(11800, new Color(0.74197871f, 0.83255582f, 1.00000000f));
            temp2color.Add(11900, new Color(0.73985682f, 0.83110912f, 1.00000000f));
            temp2color.Add(12000, new Color(0.73778012f, 0.82969211f, 1.00000000f));
            temp2color.Add(12100, new Color(0.73574723f, 0.82830393f, 1.00000000f));
            temp2color.Add(12200, new Color(0.73375683f, 0.82694373f, 1.00000000f));
            temp2color.Add(12300, new Color(0.73180765f, 0.82561071f, 1.00000000f));
            temp2color.Add(12400, new Color(0.72989845f, 0.82430410f, 1.00000000f));
            temp2color.Add(12500, new Color(0.72802807f, 0.82302316f, 1.00000000f));
            temp2color.Add(12600, new Color(0.72619537f, 0.82176715f, 1.00000000f));
            temp2color.Add(12700, new Color(0.72439927f, 0.82053539f, 1.00000000f));
            temp2color.Add(12800, new Color(0.72263872f, 0.81932722f, 1.00000000f));
            temp2color.Add(12900, new Color(0.72091270f, 0.81814197f, 1.00000000f));
            temp2color.Add(13000, new Color(0.71922025f, 0.81697905f, 1.00000000f));
            temp2color.Add(13100, new Color(0.71756043f, 0.81583783f, 1.00000000f));
            temp2color.Add(13200, new Color(0.71593234f, 0.81471775f, 1.00000000f));
            temp2color.Add(13300, new Color(0.71433510f, 0.81361825f, 1.00000000f));
            temp2color.Add(13400, new Color(0.71276788f, 0.81253878f, 1.00000000f));
            temp2color.Add(13500, new Color(0.71122987f, 0.81147883f, 1.00000000f));
            temp2color.Add(13600, new Color(0.70972029f, 0.81043789f, 1.00000000f));
            temp2color.Add(13700, new Color(0.70823838f, 0.80941546f, 1.00000000f));
            temp2color.Add(13800, new Color(0.70678342f, 0.80841109f, 1.00000000f));
            temp2color.Add(13900, new Color(0.70535469f, 0.80742432f, 1.00000000f));
            temp2color.Add(14000, new Color(0.70395153f, 0.80645469f, 1.00000000f));
            temp2color.Add(14100, new Color(0.70257327f, 0.80550180f, 1.00000000f));
            temp2color.Add(14200, new Color(0.70121928f, 0.80456522f, 1.00000000f));
            temp2color.Add(14300, new Color(0.69988894f, 0.80364455f, 1.00000000f));
            temp2color.Add(14400, new Color(0.69858167f, 0.80273941f, 1.00000000f));
            temp2color.Add(14500, new Color(0.69729688f, 0.80184943f, 1.00000000f));
            temp2color.Add(14600, new Color(0.69603402f, 0.80097423f, 1.00000000f));
            temp2color.Add(14700, new Color(0.69479255f, 0.80011347f, 1.00000000f));
            temp2color.Add(14800, new Color(0.69357196f, 0.79926681f, 1.00000000f));
            temp2color.Add(14900, new Color(0.69237173f, 0.79843391f, 1.00000000f));
            temp2color.Add(15000, new Color(0.69119138f, 0.79761446f, 1.00000000f));
            temp2color.Add(15100, new Color(0.69003044f, 0.79680814f, 1.00000000f));
            temp2color.Add(15200, new Color(0.68888844f, 0.79601466f, 1.00000000f));
            temp2color.Add(15300, new Color(0.68776494f, 0.79523371f, 1.00000000f));
            temp2color.Add(15400, new Color(0.68665951f, 0.79446502f, 1.00000000f));
            temp2color.Add(15500, new Color(0.68557173f, 0.79370830f, 1.00000000f));
            temp2color.Add(15600, new Color(0.68450119f, 0.79296330f, 1.00000000f));
            temp2color.Add(15700, new Color(0.68344751f, 0.79222975f, 1.00000000f));
            temp2color.Add(15800, new Color(0.68241029f, 0.79150740f, 1.00000000f));
            temp2color.Add(15900, new Color(0.68138918f, 0.79079600f, 1.00000000f));
            temp2color.Add(16000, new Color(0.68038380f, 0.79009531f, 1.00000000f));
            temp2color.Add(16100, new Color(0.67939381f, 0.78940511f, 1.00000000f));
            temp2color.Add(16200, new Color(0.67841888f, 0.78872517f, 1.00000000f));
            temp2color.Add(16300, new Color(0.67745866f, 0.78805526f, 1.00000000f));
            temp2color.Add(16400, new Color(0.67651284f, 0.78739518f, 1.00000000f));
            temp2color.Add(16500, new Color(0.67558112f, 0.78674472f, 1.00000000f));
            temp2color.Add(16600, new Color(0.67466317f, 0.78610368f, 1.00000000f));
            temp2color.Add(16700, new Color(0.67375872f, 0.78547186f, 1.00000000f));
            temp2color.Add(16800, new Color(0.67286748f, 0.78484907f, 1.00000000f));
            temp2color.Add(16900, new Color(0.67198916f, 0.78423512f, 1.00000000f));
            temp2color.Add(17000, new Color(0.67112350f, 0.78362984f, 1.00000000f));
            temp2color.Add(17100, new Color(0.67027024f, 0.78303305f, 1.00000000f));
            temp2color.Add(17200, new Color(0.66942911f, 0.78244457f, 1.00000000f));
            temp2color.Add(17300, new Color(0.66859988f, 0.78186425f, 1.00000000f));
            temp2color.Add(17400, new Color(0.66778228f, 0.78129191f, 1.00000000f));
            temp2color.Add(17500, new Color(0.66697610f, 0.78072740f, 1.00000000f));
            temp2color.Add(17600, new Color(0.66618110f, 0.78017057f, 1.00000000f));
            temp2color.Add(17700, new Color(0.66539706f, 0.77962127f, 1.00000000f));
            temp2color.Add(17800, new Color(0.66462376f, 0.77907934f, 1.00000000f));
            temp2color.Add(17900, new Color(0.66386098f, 0.77854465f, 1.00000000f));
            temp2color.Add(18000, new Color(0.66310852f, 0.77801705f, 1.00000000f));
            temp2color.Add(18100, new Color(0.66236618f, 0.77749642f, 1.00000000f));
            temp2color.Add(18200, new Color(0.66163375f, 0.77698261f, 1.00000000f));
            temp2color.Add(18300, new Color(0.66091106f, 0.77647551f, 1.00000000f));
            temp2color.Add(18400, new Color(0.66019791f, 0.77597498f, 1.00000000f));
            temp2color.Add(18500, new Color(0.65949412f, 0.77548090f, 1.00000000f));
            temp2color.Add(18600, new Color(0.65879952f, 0.77499315f, 1.00000000f));
            temp2color.Add(18700, new Color(0.65811392f, 0.77451161f, 1.00000000f));
            temp2color.Add(18800, new Color(0.65743716f, 0.77403618f, 1.00000000f));
            temp2color.Add(18900, new Color(0.65676908f, 0.77356673f, 1.00000000f));
            temp2color.Add(19000, new Color(0.65610952f, 0.77310316f, 1.00000000f));
            temp2color.Add(19100, new Color(0.65545831f, 0.77264537f, 1.00000000f));
            temp2color.Add(19200, new Color(0.65481530f, 0.77219324f, 1.00000000f));
            temp2color.Add(19300, new Color(0.65418036f, 0.77174669f, 1.00000000f));
            temp2color.Add(19400, new Color(0.65355332f, 0.77130560f, 1.00000000f));
            temp2color.Add(19500, new Color(0.65293404f, 0.77086988f, 1.00000000f));
            temp2color.Add(19600, new Color(0.65232240f, 0.77043944f, 1.00000000f));
            temp2color.Add(19700, new Color(0.65171824f, 0.77001419f, 1.00000000f));
            temp2color.Add(19800, new Color(0.65112144f, 0.76959404f, 1.00000000f));
            temp2color.Add(19900, new Color(0.65053187f, 0.76917889f, 1.00000000f));
            temp2color.Add(20000, new Color(0.64994941f, 0.76876866f, 1.00000000f));
            temp2color.Add(20100, new Color(0.64937392f, 0.76836326f, 1.00000000f));
            temp2color.Add(20200, new Color(0.64880528f, 0.76796263f, 1.00000000f));
            temp2color.Add(20300, new Color(0.64824339f, 0.76756666f, 1.00000000f));
            temp2color.Add(20400, new Color(0.64768812f, 0.76717529f, 1.00000000f));
            temp2color.Add(20500, new Color(0.64713935f, 0.76678844f, 1.00000000f));
            temp2color.Add(20600, new Color(0.64659699f, 0.76640603f, 1.00000000f));
            temp2color.Add(20700, new Color(0.64606092f, 0.76602798f, 1.00000000f));
            temp2color.Add(20800, new Color(0.64553103f, 0.76565424f, 1.00000000f));
            temp2color.Add(20900, new Color(0.64500722f, 0.76528472f, 1.00000000f));
            temp2color.Add(21000, new Color(0.64448939f, 0.76491935f, 1.00000000f));
            temp2color.Add(21100, new Color(0.64397745f, 0.76455808f, 1.00000000f));
            temp2color.Add(21200, new Color(0.64347129f, 0.76420082f, 1.00000000f));
            temp2color.Add(21300, new Color(0.64297081f, 0.76384753f, 1.00000000f));
            temp2color.Add(21400, new Color(0.64247594f, 0.76349813f, 1.00000000f));
            temp2color.Add(21500, new Color(0.64198657f, 0.76315256f, 1.00000000f));
            temp2color.Add(21600, new Color(0.64150261f, 0.76281076f, 1.00000000f));
            temp2color.Add(21700, new Color(0.64102399f, 0.76247267f, 1.00000000f));
            temp2color.Add(21800, new Color(0.64055061f, 0.76213824f, 1.00000000f));
            temp2color.Add(21900, new Color(0.64008239f, 0.76180740f, 1.00000000f));
            temp2color.Add(22000, new Color(0.63961926f, 0.76148010f, 1.00000000f));
            temp2color.Add(22100, new Color(0.63916112f, 0.76115628f, 1.00000000f));
            temp2color.Add(22200, new Color(0.63870790f, 0.76083590f, 1.00000000f));
            temp2color.Add(22300, new Color(0.63825953f, 0.76051890f, 1.00000000f));
            temp2color.Add(22400, new Color(0.63781592f, 0.76020522f, 1.00000000f));
            temp2color.Add(22500, new Color(0.63737701f, 0.75989482f, 1.00000000f));
            temp2color.Add(22600, new Color(0.63694273f, 0.75958764f, 1.00000000f));
            temp2color.Add(22700, new Color(0.63651299f, 0.75928365f, 1.00000000f));
            temp2color.Add(22800, new Color(0.63608774f, 0.75898278f, 1.00000000f));
            temp2color.Add(22900, new Color(0.63566691f, 0.75868499f, 1.00000000f));
            temp2color.Add(23000, new Color(0.63525042f, 0.75839025f, 1.00000000f));
            temp2color.Add(23100, new Color(0.63483822f, 0.75809849f, 1.00000000f));
            temp2color.Add(23200, new Color(0.63443023f, 0.75780969f, 1.00000000f));
            temp2color.Add(23300, new Color(0.63402641f, 0.75752379f, 1.00000000f));
            temp2color.Add(23400, new Color(0.63362667f, 0.75724075f, 1.00000000f));
            temp2color.Add(23500, new Color(0.63323097f, 0.75696053f, 1.00000000f));
            temp2color.Add(23600, new Color(0.63283925f, 0.75668310f, 1.00000000f));
            temp2color.Add(23700, new Color(0.63245144f, 0.75640840f, 1.00000000f));
            temp2color.Add(23800, new Color(0.63206749f, 0.75613641f, 1.00000000f));
            temp2color.Add(23900, new Color(0.63168735f, 0.75586707f, 1.00000000f));
            temp2color.Add(24000, new Color(0.63131096f, 0.75560036f, 1.00000000f));
            temp2color.Add(24100, new Color(0.63093826f, 0.75533624f, 1.00000000f));
            temp2color.Add(24200, new Color(0.63056920f, 0.75507467f, 1.00000000f));
            temp2color.Add(24300, new Color(0.63020374f, 0.75481562f, 1.00000000f));
            temp2color.Add(24400, new Color(0.62984181f, 0.75455904f, 1.00000000f));
            temp2color.Add(24500, new Color(0.62948337f, 0.75430491f, 1.00000000f));
            temp2color.Add(24600, new Color(0.62912838f, 0.75405319f, 1.00000000f));
            temp2color.Add(24700, new Color(0.62877678f, 0.75380385f, 1.00000000f));
            temp2color.Add(24800, new Color(0.62842852f, 0.75355685f, 1.00000000f));
            temp2color.Add(24900, new Color(0.62808356f, 0.75331217f, 1.00000000f));
            temp2color.Add(25000, new Color(0.62774186f, 0.75306977f, 1.00000000f));
            temp2color.Add(25100, new Color(0.62740336f, 0.75282962f, 1.00000000f));
        }
    }
}