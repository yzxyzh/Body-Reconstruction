#pragma once

#include <vtkImageActor.h>
#include <vtkContourFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkImageData.h>
#include <vtkMarchingCubes.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkPlaneSource.h>

#include <vtkCylinderSource.h>

#include <opencv2/core/core.hpp>
#include "myInteractorStyle.h"

#include <vtkLineSource.h>
#include <vtkRegularPolygonSource.h>

#include <vtkAxisActor.h>

#include <vtkOrientationMarkerWidget.h>
#include <vtkAnnotatedCubeActor.h>
#include<vtkDistanceWidget.h>
#include<vtkDistanceRepresentation3D.h>

#include<vtkUnsignedCharArray.h>
//#include <vtkTextActor.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

using namespace cv;

class Pool{
public:
    static vtkSmartPointer<vtkRenderer>                       renderer;
    static vtkSmartPointer<vtkRenderWindow>                   renderWindow;
    static vtkSmartPointer<vtkRenderWindowInteractor>         interactor;
    static vtkSmartPointer<vtkCamera>                         camera;
    static vtkSmartPointer<vtkSphereSource>                   sphere;
    static vtkSmartPointer<vtkActor>                          sphereActor;
    static vtkSmartPointer<vtkPolyDataMapper>                 sphereMapper;
    static vtkSmartPointer<vtkSphereSource>                   Point;
    static vtkSmartPointer<vtkActor>                          PointActor;
    static vtkSmartPointer<vtkPolyDataMapper>                 PointMapper;
    static vtkSmartPointer<vtkTransform>                      userTransform;
    static vtkSmartPointer<vtkCylinderSource>                 cylinder;
    static vtkSmartPointer<vtkPolyDataMapper>                 cylinderMapper;
    static vtkSmartPointer<vtkActor>                          cylinderActor;
    static vtkSmartPointer<vtkLineSource>                     line;
    static vtkSmartPointer<vtkPolyDataMapper>                 lineMapper;
    static vtkSmartPointer<vtkActor>                          lineActor;
    static vtkMatrix4x4*                                      cameraTransform;
    
    static myInteractorStyle*                                 trackball;
    static Point3d                                            currPoint;
    static Point3d                                            sphereCenter;
    static Point3d                                            tumorCenter;
    static vector<vector<Point3d> >                           interSectors;
    static vector<vector<double> >                            colors;
    static vector<vtkSmartPointer<vtkActor> >                 lastIntersectionActors;
	//traverse Point2 made by CRD
	static int                                                 ballNumber;
	static vector<Point3d>                                     tmpCenter;
    
    static vector<vtkSmartPointer<vtkRegularPolygonSource> >  bigSpheres;
    static vector<vtkSmartPointer<vtkPolyDataMapper> >        bigSpheresMapper;
    static vector<vtkSmartPointer<vtkActor> >                 bigSpheresActor;
    
    //Scene variables
    static vector<vtkSmartPointer<vtkActor> >                 sceneActors;
    static vector<vtkSmartPointer<vtkPolyData> >              scenePolyDatas;
    
    static vector<vtkSmartPointer<vtkAxisActor> >             axisActors;
    
	//add axis
	static vtkSmartPointer<vtkAnnotatedCubeActor>             mainAxisActor;
	static vtkSmartPointer<vtkOrientationMarkerWidget>        mainAxisWidget;
	//add distance measure
	static vtkSmartPointer<vtkDistanceWidget>                 distanceWidget;
	static vtkSmartPointer<vtkDistanceRepresentation3D>       rep3D;
	     
              
	//multiball
	static vector<vtkSmartPointer<vtkSphereSource> >          coverBallSources;
	static vector<vtkSmartPointer<vtkPolyDataMapper> >        coverBallMappers;
	static vector<vtkSmartPointer<vtkActor> >                 coverBallActors;

	static vector<vtkSmartPointer<vtkSphereSource> >          coverCenterSources;
	static vector<vtkSmartPointer<vtkPolyDataMapper> >        coverCenterMappers;
	static vector<vtkSmartPointer<vtkActor> >                 coverCenterActors;

	

	static vtkSmartPointer<vtkTextActor>                      thetaActor;
	static vtkSmartPointer<vtkTextActor>                      phiActor;
	
	static int                                                currentBallIndex;

    static void                                               RefreshScene();
};










