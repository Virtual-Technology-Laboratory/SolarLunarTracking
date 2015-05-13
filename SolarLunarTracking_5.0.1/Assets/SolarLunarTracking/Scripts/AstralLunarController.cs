/*
 * Copyright (c) 2014, Roger Lew (rogerlew.gmail.com)
 * Date: 2/25/2015 - 5/13/2015
 * License: BSD (3-clause license)
 * 
 * The project described was supported by NSF award number IIA-1301792
 * from the NSF Idaho EPSCoR Program and by the National Science Foundation.
 * 
 */


using System;
using UnityEngine;
using System.Collections;
using System.Collections.Generic;

using VTL.SimTimeControls;

namespace VTL.SolarLunarTracking
{
    public class AstralLunarController : MonoBehaviour
    {
        private float latitude = 43f;
        private float longitude = -116f;
        private float altitude = 1160;
        private float timeZone = -7;

//        private Light moon;
        private GameObject plane;
        private List<Texture> moonTextures;
        private int currentTex = 11;

        private GameObject slider;
        private TimeSlider TimeSliderScript;

        private Astral astral;
        private Lunar lunar;
        public Lunar lunarProperty
        {
            get { return lunar; }
            set { lunar = value; }
        }

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
            TimeSliderScript = (TimeSlider)slider.GetComponent(typeof(TimeSlider));

            astral = new Astral(longitude, latitude, altitude, timeZone);
//            moon = transform.Find("Plane").gameObject.GetComponent<Light>();
            plane = transform.Find("Plane").gameObject;

            moonTextures = new List<Texture>();
            for (int i = 0; i < 28; i++)
            {
                string fname = string.Format("Assets/SolarLunarTracking/moonTextures/Moon{0}.png", i.ToString("D2"));
                moonTextures.Add(Resources.LoadAssetAtPath<Texture>(fname));
            }
        }



        void Update()
        {
            DateTime datetime = TimeSliderScript.SimTime;

            lunar = astral.lunar(datetime);
            transform.eulerAngles = new Vector3(270 + (float)lunar.elevation, 180 + (float)lunar.azimuth, 0);

            int matInd = Mathf.FloorToInt((float)lunar.phase * 28f);

            if (lunar.elevation < -20)
            {
                if (currentTex != matInd)
                {
                    plane.GetComponent<Renderer>().material.mainTexture = moonTextures[matInd];
                    currentTex = matInd;
                }
            }

        }
    }
}