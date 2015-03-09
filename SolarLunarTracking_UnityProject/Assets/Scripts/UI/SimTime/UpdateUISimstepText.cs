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

public class UpdateUISimstepText : MonoBehaviour {

    private GameObject slider;
    private TimeSlider timeSlider;

    private Text text;

    // Use this for initialization
    void Start()
    {
        slider = GameObject.FindWithTag("TimeSlider");
        timeSlider = (TimeSlider)slider.GetComponent(typeof(TimeSlider));
        
        text = GetComponent<Text>();
    }
	
	// Update is called once per frame
	void Update () {
        text.text = timeSlider.SimStep.ToString("D4");
	}
}
