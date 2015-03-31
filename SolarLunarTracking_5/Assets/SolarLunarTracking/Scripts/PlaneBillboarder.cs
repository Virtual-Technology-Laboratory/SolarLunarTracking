using UnityEngine;
using System.Collections;

public class PlaneBillboarder : MonoBehaviour {


    void Start()
    {
        GetComponent<Renderer>().material.renderQueue = 3000;
    }

	// Update is called once per frame
    void Update()
    {
        transform.LookAt(Camera.main.transform.position, Vector3.up);
        transform.Rotate(new Vector3(90, 0, 0));
	}
}
