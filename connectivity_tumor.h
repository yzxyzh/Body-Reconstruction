#pragma once
#include<vector>
#include<iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
#ifndef _in_
#define _in_
#endif 

#ifndef _out_
#define _out_
#endif



	namespace connectivity{

		
		void ComputeConnectedParts(_in_ Mat * Mask,_in_ int dims1,_in_ int dims2,_in_ int dims3,_out_ vector<vector<Point3i> >& connectedMap );

	}
