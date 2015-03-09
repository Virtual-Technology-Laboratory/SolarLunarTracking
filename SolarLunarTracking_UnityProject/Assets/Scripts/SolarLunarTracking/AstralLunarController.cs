using System;
using UnityEngine;
using System.Collections;


public class AstralLunarController : MonoBehaviour 
{
    public float latitude = 43f;
    public float longitude = -116f;
    public float altitude = 1160;
    public float timeZone = -7;

    private Light moon;

    private GameObject slider;
    private TimeSlider TimeSliderScript;

    private Astral astral;
    private Lunar lunar;
    public Lunar lunarProperty
    {
        get { return lunar; }
        set { lunar = value; }
    }

    // Use this for initialization
    void Start()
    {
        slider = GameObject.FindWithTag("TimeSlider");
        TimeSliderScript = (TimeSlider)slider.GetComponent(typeof(TimeSlider));

        astral = new Astral(longitude, latitude, altitude, timeZone);
        foreach (Transform t in transform)
         {
             if (t.name == "Light")
                 moon = t.GetComponent<Light>();
         }
    }

    // Update is called once per frame
    void Update()
    {
        DateTime datetime = TimeSliderScript.SimTime;
        
        lunar = astral.lunar(datetime);

//        Debug.Log(String.Format("Lunar: {0}, {1}", lunar_az, lunar_el));

        transform.eulerAngles = new Vector3(270 + (float)lunar.elevation, 180 + (float)lunar.azimuth, 0);

    }
}