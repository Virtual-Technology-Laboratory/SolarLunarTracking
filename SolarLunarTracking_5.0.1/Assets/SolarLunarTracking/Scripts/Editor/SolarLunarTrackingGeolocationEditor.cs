using UnityEngine;
using UnityEditor;
using System.Collections;

namespace VTL.SolarLunarTracking
{
    [CustomEditor(typeof(SolarLunarTrackingGeolocation))]
    public class SolarLunarTrackingGeolocationEditor : Editor
    {
        public override void OnInspectorGUI()
        {
            DrawDefaultInspector();
            SolarLunarTrackingGeolocation script = (SolarLunarTrackingGeolocation)target;
            if (GUILayout.Button("Drop on Terrain"))
                script.DropOnTerrain();
        }
    }
}