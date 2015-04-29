#pragma once
#include"vtkMetaImageReader.h"
#include"vtkImageData.h"
#include"vtkImageActor.h"
#include"vtkRenderer.h"
#include"vtkRenderWindow.h"
#include"vtkRenderWindowInteractor.h"
#include"vtkInteractorStyleImage.h"
#include"vtkActor.h"
#include"vtkMarchingCubes.h"
#include"vtkVectorNorm.h"
#include"vtkDataSetMapper.h"
#include"vtkProperty.h"
#include"vtkSmoothPolyDataFilter.h"
#include"vtkBoxWidget.h"
#include"vtkPolyDataMapper.h"
#include"vtkPolyDataNormals.h"
#include"vtkPLYWriter.h"
#include"vtkPLYReader.h"
#include"vtkQuadricClustering.h"
#include"vtkRegularPolygonSource.h"
#include"vtkCamera.h"
#include<iostream>
#include<string>
#include<math.h>
using namespace std;
#ifndef _In_
#define _In_
#endif // !_In_

#ifndef _Out_
#define _Out_
#endif

namespace CRD{
	void drawer(_In_ const string& skin_mhd_file_name);
}