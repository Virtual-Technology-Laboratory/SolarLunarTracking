/*
 * Copyright (c) 2014, Roger Lew (rogerlew.gmail.com)
 * Date: 2/5/2015
 * License: BSD (3-clause license)
 * 
 * The project described was supported by NSF award number IIA-1301792
 * from the NSF Idaho EPSCoR Program and by the National Science Foundation.
 * 
 */

using UnityEngine;
using System.Collections;
using UnityEngine.UI;

public class PlaySpeedSliderControl : MonoBehaviour {

    private GameObject slider;
    private TimeSlider timeSlider;

    private Slider speedControlSlider;

    private Text speedText;

	// Use this for initialization
	void Start () 
    {
        slider = GameObject.FindWithTag("TimeSlider");
        timeSlider = (TimeSlider)slider.GetComponent(typeof(TimeSlider));

        speedControlSlider = GameObject.FindWithTag("SpeedControlSlider").GetComponent<Slider>();
        speedText = GameObject.FindWithTag("SpeedText").GetComponent<Text>();

        OnValueChange(3);
	}
	
	// Update is called once per frame
	void Update () 
    {
	}
    
    public void OnValueChange()
    {
        OnValueChange(speedControlSlider.value);
    }

    public void OnValueChange(float value)
    {
        timeSlider.PlaySpeed = Mathf.Pow(2f, value) / 8f;
        float spd = timeSlider.PlaySpeed;
        if (spd >= 1)
            speedText.text = string.Format("1s:{0}h", (int)spd);
        else
            speedText.text = string.Format("1s:{0}m", (float)(spd*60f));

    }

}
