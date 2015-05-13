/*
 * Copyright (c) 2014, Roger Lew (rogerlew.gmail.com)
 * Date: 2/25/2015 - 5/13/2015
 * License: BSD (3-clause license)
 * 
 * The project described was supported by NSF award number IIA-1301792
 * from the NSF Idaho EPSCoR Program and by the National Science Foundation.
 * 
 */

#define NODEBUG

using System;
using UnityEngine;
using System.Collections;

using VTL.SimTimeControls;

namespace VTL.SolarLunarTracking
{
    public class AstralSunController : MonoBehaviour
    {

        private float latitude = 43f;
        private float longitude = -116f;
        private float altitude = 1160;
        private float timeZone = -7;

        private Light sun;

        private GameObject slider;
        private TimeSlider timeSlider;

        private Astral astral;
        private Solar solar;
        public Solar solarProperty
        {
            get { return solar; }
            set { solar = value; }
        }


        // Use this for initialization
        void Start()
        {
            // Read the geolocation from the parent
            GameObject parent = transform.parent.gameObject;
            SolarLunarTrackingGeolocation geoLoc = parent.GetComponent<SolarLunarTrackingGeolocation>();
            latitude = geoLoc.latitude;
            longitude = geoLoc.longitude;
            altitude = geoLoc.altitude;
            timeZone = geoLoc.timeZone;

#if DEBUG
            Debug.Log(String.Format("Location: {0}, {1}", longitude, latitude));
            Debug.Log(String.Format("Altitude: {0}", altitude));
            Debug.Log(String.Format("Time Zone: {0}", timeZone));
#endif

            slider = GameObject.FindWithTag("TimeSlider");
            timeSlider = (TimeSlider)slider.GetComponent<TimeSlider>();

            astral = new Astral(longitude, latitude, altitude, timeZone);
            foreach (Transform t in transform)
            {
                if (t.name == "Light")
                    sun = t.GetComponent<Light>();
            }
        }

        // Update is called once per frame
        void Update()
        {

            DateTime datetime = timeSlider.SimTime;
            solar = astral.solar(datetime);

#if DEBUG
            Debug.Log(datetime);
            Debug.Log(String.Format("Solar: {0}, {1}", solar.azimuth, solar.elevation));
#endif
            transform.eulerAngles = new Vector3(270 + (float)solar.elevation, 180 + (float)solar.azimuth, 0);
            sun.color = solar.color_temp;

        }
    }
}