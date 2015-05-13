/*
 * Copyright (c) 2014, Roger Lew (rogerlew.gmail.com)
 * Date: 2/25/2015 - 5/13/2015
 * License: BSD (3-clause license)
 * 
 * The project described was supported by NSF award number IIA-1301792
 * from the NSF Idaho EPSCoR Program and by the National Science Foundation.
 * 
 */

using UnityEngine;
using System.Collections;

namespace VTL.SolarLunarTracking
{
    public class PlaneBillboarder : MonoBehaviour
    {

        void Start()
        {
            GetComponent<Renderer>().material.renderQueue = 3000;
        }

        void Update()
        {
            transform.LookAt(Camera.main.transform.position, Vector3.up);
            transform.Rotate(new Vector3(90, 0, 0));
        }
    }
}