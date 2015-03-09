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

public class TimeCtrlToggleGroupEventHandlers : MonoBehaviour 
{
    ToggleGroup toggleGroup;

    private GameObject slider;
    private TimeSlider timeSlider;

	// Use this for initialization
    void Start()
    {
        slider = GameObject.FindWithTag("TimeSlider"); 
        timeSlider = (TimeSlider)slider.GetComponent(typeof(TimeSlider));
	}
	
	// Called by UI Event Handler
	public void onValueChange(Toggle playToggle) {
        timeSlider.IsPlaying = playToggle.isOn;
	}
}
