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
    public class SolarLunarTrackingGeolocation : MonoBehaviour
    {

        public float latitude = 43f;
        public float longitude = -116f;
        public float altitude = 1160;
        public float timeZone = -7;

        public void DropOnTerrain()
        {
           RaycastHit hitInfo;
            if (Physics.Raycast(transform.position, Vector3.down, out hitInfo, 1000000))
            {
                var newpos = transform.position;
                newpos.y -= hitInfo.distance;
                transform.position = newpos;
            }
            else
            {
                Debug.Log("Raycast did not hit Terrain");
            }
        }
    }
}
