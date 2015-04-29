#pragma once
#include"vtkSphereSource.h"
#include"vtkPolyDataMapper.h"
#include"vtkActor.h"
#include"vtkRenderer.h"
#include"vtkRenderWindow.h"
#include"vtkRenderWindowInteractor.h"
#include"vtkMetaImageReader.h"
#include"vtkImageActor.h"
#include"vtkImageData.h"
#include"vtkProperty.h"
#include"vtkPoints.h"
#include"vtkCellArray.h"
#include"vtkPolyData.h"
#include <vector>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <string>

using namespace std;
using namespace cv;
#ifndef _In_
#define _In_
#endif // !_In_

#ifndef _Out_
#define _Out_
#endif



namespace YZX{
	typedef vector<vector<double> > CenterType;
	typedef vector<double> RadiusType;

	void balls_Execute(
		_In_ const string& tumor_mhd_file_name,//进入的肿瘤文件名称
		_In_ const double needle_radius,//针半径
		_Out_ int& min_balls_num,//输出的最小覆盖球个数
		_Out_ CenterType& centers,//输出的球中心
		_Out_ RadiusType& radiuses//输出的球半径
		);

	void GetBallCover(
		_In_ const Mat& PointSeries,
		_In_ const int clusterNum,
		_Out_ CenterType& out_centers,
		_Out_ RadiusType& out_radiuses,
		_Out_ double& maxR
		);
}