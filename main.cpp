#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

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

#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"

#include <vtkProperty.h>
#include <vtkStringArray.h>

#include "myInteractorStyle.h"

#include "Pool.h"

#include <vtkOBBTree.h>

#include <random>
#include "drawing.h"

#include "Cover.h"


using namespace cv;

void AppendScene(const string& inFile,
                 const double opacity,
                 const double r,
                 const double g,
                 const double b,
                 bool  isTumor = false)
{
    vtkSmartPointer< vtkMetaImageReader> reader = vtkSmartPointer<vtkMetaImageReader>::New();
    reader->SetFileName(inFile.c_str());
    reader->Update();
    
    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
    image=reader->GetOutput();
    
    //image->Print(cout);
    
    vtkSmartPointer<vtkMarchingCubes> iso = vtkSmartPointer<vtkMarchingCubes>::New();
#if VTK_MAJOR_VERSION <6
    iso->SetInput(reader->GetOutput());
#else
    iso->SetInputData(reader->GetOutput());
#endif
    iso->SetNumberOfContours(1);
    iso->SetValue(0,1);
    iso->ComputeGradientsOn();
    iso->ComputeNormalsOn();
    iso->ComputeScalarsOff();
    
    
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter =
    vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputConnection(iso->GetOutputPort());
    smoothFilter->SetNumberOfIterations(10);
    smoothFilter->FeatureEdgeSmoothingOn();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();
    
    vtkSmartPointer<vtkQuadricClustering> decimate = vtkSmartPointer<vtkQuadricClustering>::New();
    decimate->SetNumberOfXDivisions(1200);
    decimate->SetNumberOfYDivisions(1200);
    decimate->SetNumberOfZDivisions(1200);
    
#if VTK_MAJOR_VERSION <6
    decimate->SetInput(smoothFilter->GetOutput());
#else
    decimate->SetInputData(smoothFilter->GetOutput());
#endif
    decimate->Update();
    
    vtkPolyData* pd = decimate->GetOutput();
    Pool::scenePolyDatas.push_back(pd);
    
    vtkSmartPointer<vtkPolyDataMapper> dataMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    dataMapper->SetInputConnection(decimate->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(dataMapper);
    actor->GetProperty()->SetOpacity(opacity);
    actor->GetProperty()->SetColor(r, g, b);

	



    Pool::sceneActors.push_back(actor);
    
    
    
    Pool::trackball->AddObject(pd);
    
    
    if(isTumor)
    {

        double* currBound = pd->GetBounds();
            
	    Pool::tumorCenter.x = 0.5*(currBound[0]+currBound[1]);
		Pool::tumorCenter.y = 0.5*(currBound[2]+currBound[3]);
		Pool::tumorCenter.z = 0.5*(currBound[4]+currBound[5]);

		int ballNum = 0;
		YZX::CenterType centers;
		YZX::RadiusType radiuses;

		YZX::balls_Execute(inFile, 15.0, ballNum, centers, radiuses);

		cout<<"we have calculated : "<<ballNum<<" balls"<<endl;
		
		//travesal Point2 made by CRD
		Pool::ballNumber=ballNum;
		Pool::tmpCenter.resize(Pool::ballNumber);

		Pool::tmpCenter[0].x=0.5*(currBound[0]+currBound[1]);
		Pool::tmpCenter[0].y=0.5*(currBound[2]+currBound[3]);
		Pool::tmpCenter[0].z=0.5*(currBound[4]+currBound[5]);


		//Build new balls;
		Pool::coverBallSources.resize(ballNum);
		Pool::coverBallMappers.resize(ballNum);
		Pool::coverBallActors.resize(ballNum);

		Pool::coverCenterSources.resize(ballNum);
		Pool::coverCenterMappers.resize(ballNum);
		Pool::coverCenterActors.resize(ballNum);
		for (int i=0; i<ballNum; i++) {
			Pool::coverBallSources[i] = vtkSmartPointer<vtkSphereSource>::New();
			double xC = centers[i][0];
			double yC = centers[i][1];
			double zC = centers[i][2];
			double rC = radiuses[i];

			cout<<i<<" th ball center = "<<xC<<","<<yC<<","<<zC<<" radius = "<<rC<<endl;
			//made by CRD
			Pool::tmpCenter[i].x=xC;
			Pool::tmpCenter[i].y=yC;
			Pool::tmpCenter[i].z=zC;

			Pool::coverBallSources[i]->SetCenter(xC, yC, zC);
			Pool::coverBallSources[i]->SetRadius(15.0);   //this is the radius of the injection
			Pool::coverBallSources[i]->SetThetaResolution(30);
			Pool::coverBallSources[i]->SetPhiResolution(30);
			Pool::coverBallSources[i]->Update();

			Pool::coverBallMappers[i] = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION < 6
			Pool::coverBallMappers[i]->SetInput(Pool::coverBallSources[i]->GetOutput());
#else
            Pool::coverBallMappers[i]->SetInputData(Pool::coverBallSources[i]->GetOutput());
#endif
			Pool::coverBallActors[i] = vtkSmartPointer<vtkActor>::New();
			Pool::coverBallActors[i]->SetMapper(Pool::coverBallMappers[i]);
			Pool::coverBallActors[i]->GetProperty()->SetColor(0.0, 1.0, 1.0);
			Pool::coverBallActors[i]->GetProperty()->SetOpacity(0.95);

			Pool::coverCenterSources[i] = vtkSmartPointer<vtkSphereSource>::New();
			Pool::coverCenterSources[i]->SetCenter(xC, yC, zC);
			Pool::coverCenterSources[i]->SetRadius(2.0);
			/*Pool::coverCenterSources[i]->SetThetaResolution(30);
			Pool::coverCenterSources[i]->SetPhiResolution(30);
			Pool::coverCenterSources[i]->Update();*/

			Pool::coverCenterMappers[i] = vtkSmartPointer<vtkPolyDataMapper>::New();
            
#if VTK_MAJOR_VERSION < 6
			Pool::coverCenterMappers[i]->SetInput(Pool::coverCenterSources[i]->GetOutput());
#else
            Pool::coverCenterMappers[i]->SetInputData(Pool::coverCenterSources[i]->GetOutput());
#endif
			

			Pool::coverCenterActors[i] = vtkSmartPointer<vtkActor>::New();
			Pool::coverCenterActors[i]->SetMapper(Pool::coverCenterMappers[i]);
			Pool::coverCenterActors[i]->GetProperty()->SetColor(1.0, 1.0, 0.0);
			Pool::coverCenterActors[i]->GetProperty()->SetOpacity(1.0);
		}

		Pool::tumorCenter.x = centers[0][0];
		Pool::tumorCenter.y = centers[0][1];
		Pool::tumorCenter.z = centers[0][2];
		

    }
 
}


void SetupWindow()
{

    Pool::userTransform->Identity();
    
    Pool::renderer ->SetBackground(0.0, 0.0, 0.0);
    
    Pool::renderWindow->AddRenderer(Pool::renderer);
    Pool::renderWindow->SetSize(600, 600);
    
    Pool::interactor->SetRenderWindow(Pool::renderWindow);
//    Pool::cylinder->SetCenter(Pool::tumorCenter.x, Pool::tumorCenter.y, Pool::tumorCenter.z);
//    Pool::cylinder->SetRadius(0.3);
//    Pool::cylinder->SetHeight(0.5);
//    Pool::cylinder->SetResolution(100);
//    Pool::cylinder->Update();
//    
//    
//    Pool::cylinderMapper->SetInputConnection(Pool::cylinder->GetOutputPort());
//    //Pool::cylinderMapper->SetInputData(Pool::cylinder->GetOutput());
//    Pool::cylinderActor->SetMapper(Pool::cylinderMapper);
//    Pool::cylinderActor->GetProperty()->SetOpacity(0.7);
    

    
    

}

void SetupAABB()
{
    double xMin = 9999999;
    double xMax = -9999999;
    double yMin = 9999999;
    double yMax = -9999999;
    double zMin = 9999999;
    double zMax = -9999999;
    for (int i=0; i<Pool::scenePolyDatas.size(); i++) {
        
        double* currBound = Pool::scenePolyDatas[i]->GetBounds();
        
        if(currBound[0]<xMin) xMin = currBound[0];
        if(currBound[2]<yMin) yMin = currBound[2];
        if(currBound[4]<zMin) zMin = currBound[4];
        
        if(currBound[1]>xMax) xMax = currBound[1];
        if(currBound[3]>yMax) yMax = currBound[3];
        if(currBound[5]>zMax) zMax = currBound[5];
        
    }
    
    //cout<<"AABB = ["<<xMin<<","<<xMax<<"]x["<<yMin<<","<<yMax<<"]x["<<zMin<<","<<zMax<<"]"<<endl;
    
    double radius;
    double xcenter,ycenter,zcenter;
    
    //FIXME 这边有问题，需要让他等于_xcenter,_ycenter,_zcenter
    xcenter=Pool::tumorCenter.x;
    ycenter=Pool::tumorCenter.y;
    zcenter=Pool::tumorCenter.z;
    
    radius=0.5*sqrt((xMax-xMin)*(xMax-xMin)+(yMax-yMin)*(yMax-yMin)+(zMax-zMin)*(zMax-zMin));
    
   // cout<<"center = "<<xcenter<<" "<<ycenter<<" "<<zcenter<<endl;
    //cout<<"radius = "<<radius<<endl;
    
    Pool::sphereCenter.x = xcenter;
    Pool::sphereCenter.y = ycenter;
    Pool::sphereCenter.z = zcenter;
    
    Pool::currPoint.x = xcenter+radius*1.5;
    Pool::currPoint.y = ycenter;
    Pool::currPoint.z = zcenter;
	
    
    Pool::tumorCenter.x = xcenter;
    Pool::tumorCenter.y = ycenter;
    Pool::tumorCenter.z = zcenter;
    /********INITIALIZING CAMERA*****************/
    Pool::camera->SetViewUp (0, 1, 0);
    //Pool::camera->SetPosition (0, 0, radius*3);
    //Pool::camera->SetDistance(radius*10);
    Pool::camera->SetFocalPoint (Pool::sphereCenter.x, Pool::sphereCenter.y, Pool::sphereCenter.z);
    Pool::camera->ParallelProjectionOn();
    
    /********INITIALIZING BIG SPHERE*************/
    Pool::sphere->SetCenter(Pool::sphereCenter.x, Pool::sphereCenter.y, Pool::sphereCenter.z);
    Pool::sphere->SetRadius(radius*0.8);
    Pool::sphere->SetThetaResolution(50);
    Pool::sphere->SetPhiResolution(50);
    
    Pool::sphereMapper->SetInputConnection(Pool::sphere->GetOutputPort());
    Pool::sphereActor->SetMapper(Pool::sphereMapper);
    Pool::sphereActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    Pool::sphereActor->GetProperty()->SetOpacity(0.0);
    
    
    /********INITIALIZING MOVING POINT************/
    Pool::Point->SetCenter(Pool::currPoint.x, Pool::currPoint.y, Pool::currPoint.z);
    Pool::Point->SetRadius(2.0);
    
    Pool::PointMapper->SetInputConnection(Pool::Point->GetOutputPort());
    Pool::PointActor->SetMapper(Pool::PointMapper);
    Pool::PointActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
    
    /********INITIALIZING SEG LINE*****************/
    Pool::line->SetPoint1(Pool::currPoint.x, Pool::currPoint.y, Pool::currPoint.z);
    Pool::line->SetPoint2(Pool::tumorCenter.x, Pool::tumorCenter.y, Pool::tumorCenter.z);
    Pool::lineMapper->SetInputConnection(Pool::line->GetOutputPort());
    Pool::lineActor->SetMapper(Pool::lineMapper);
    Pool::lineActor->GetProperty()->SetColor(0.0, 0.0, 1.0);
    Pool::lineActor->GetProperty()->SetLineWidth(3.0);
    
    /********INITIALIZING 3 BIG SPHERE FRAMES**************/
    Pool::bigSpheres[0] = vtkSmartPointer<vtkRegularPolygonSource>::New();
    Pool::bigSpheres[0]->GeneratePolygonOff();
    Pool::bigSpheres[0]->SetNumberOfSides(50);
    Pool::bigSpheres[0]->SetRadius(radius*1.5);
    Pool::bigSpheres[0]->SetCenter(xcenter,ycenter,zcenter);
    Pool::bigSpheres[0]->SetNormal(0,0,1);
    Pool::bigSpheres[0]->Update();   
    
    Pool::bigSpheres[1] = vtkSmartPointer<vtkRegularPolygonSource>::New();
    Pool::bigSpheres[1]->GeneratePolygonOff();
    Pool::bigSpheres[1]->SetNumberOfSides(50);
    Pool::bigSpheres[1]->SetRadius(radius*1.5);
    Pool::bigSpheres[1]->SetCenter(xcenter,ycenter,zcenter);
    Pool::bigSpheres[1]->SetNormal(0,1,0);
    Pool::bigSpheres[1]->Update();
    
    Pool::bigSpheres[2] = vtkSmartPointer<vtkRegularPolygonSource>::New();
    Pool::bigSpheres[2]->GeneratePolygonOff();
    Pool::bigSpheres[2]->SetNumberOfSides(50);
    Pool::bigSpheres[2]->SetRadius(radius*1.5);
    Pool::bigSpheres[2]->SetCenter(xcenter,ycenter,zcenter);
    Pool::bigSpheres[2]->SetNormal(1,0,0);
    Pool::bigSpheres[2]->Update();
    
    for (int i=0; i<3; i++) {
        Pool::bigSpheresMapper[i] = vtkSmartPointer<vtkPolyDataMapper>::New();
        Pool::bigSpheresActor[i] = vtkSmartPointer<vtkActor>::New();
        Pool::bigSpheresMapper[i]->SetInputConnection(Pool::bigSpheres[i]->GetOutputPort());
        Pool::bigSpheresActor[i]->SetMapper(Pool::bigSpheresMapper[i]);
    }
    
    /*********BUILD X-Y-Z AXIS***********/
    for (int i=0; i<3; i++) {
        Pool::axisActors[i] = vtkSmartPointer<vtkAxisActor>::New();
        Pool::axisActors[i]->SetLabelScale(0.2);
    }
    
    Pool::axisActors[0]->SetPoint1(xcenter,ycenter,zcenter);
    Pool::axisActors[0]->SetPoint2(xcenter+radius, ycenter, zcenter);
    vtkSmartPointer<vtkStringArray> labelsX = vtkSmartPointer<vtkStringArray>::New();
       labelsX->SetNumberOfTuples(1);
       labelsX->SetValue(0,"X");
    Pool::axisActors[0]->SetLabels(labelsX);
    
    Pool::axisActors[1]->SetPoint1(xcenter,ycenter,zcenter);
    Pool::axisActors[1]->SetPoint2(xcenter,ycenter+radius,zcenter);
    vtkSmartPointer<vtkStringArray> labelsY = vtkSmartPointer<vtkStringArray>::New();
    labelsY->SetNumberOfTuples(1);
    labelsY->SetValue(0,"Y");
    Pool::axisActors[1]->SetLabels(labelsY);
    
    Pool::axisActors[2]->SetPoint1(xcenter,ycenter,zcenter);
    Pool::axisActors[2]->SetPoint2(xcenter,ycenter,zcenter+radius);
    vtkSmartPointer<vtkStringArray> labelsZ = vtkSmartPointer<vtkStringArray>::New();
    labelsZ->SetNumberOfTuples(1);
    labelsZ->SetValue(0,"Z");
    Pool::axisActors[2]->SetLabels(labelsZ);
    
    
    
    /*********INITIALIZING INTERSECTION COLORS**********/
    std::mt19937 mtGen(0);
    std::uniform_real_distribution<double> uni(0.0,1.0);
    Pool::colors.resize(100);
    for (int i=0; i<100; i++) {
        Pool::colors[i].resize(3);
        for (int j=0; j<3; j++) {
            Pool::colors[i][j] = uni(mtGen);
        }
    }
    /*********INITIALIZING RENDERER****************/
    Pool::renderer->AddActor(Pool::lineActor);
    Pool::renderer->SetActiveCamera(Pool::camera);
    Pool::renderer->ResetCamera();
    //Pool::renderer->GetActiveCamera()->SetPosition(0, 0, radius*3);
    Pool::trackball->SetsphereActor(Pool::sphereActor);
    Pool::interactor->SetInteractorStyle(Pool::trackball);
    //Pool::renderer->AddActor(Pool::sphereActor);
    Pool::renderer->AddActor(Pool::PointActor);
    
    /****************INITIALIZING MAIN AXIS**************/
	Pool::mainAxisActor->SetXPlusFaceText("A");
	Pool::mainAxisActor->SetXMinusFaceText("P");
	Pool::mainAxisActor->SetYPlusFaceText("L");
	Pool::mainAxisActor->SetYMinusFaceText("R");
	Pool::mainAxisActor->SetZPlusFaceText("S");
	Pool::mainAxisActor->SetZMinusFaceText("I");

	Pool::mainAxisWidget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
	Pool::mainAxisWidget->SetOrientationMarker( Pool::mainAxisActor );
	Pool::mainAxisWidget->SetInteractor( Pool::interactor );
	Pool::mainAxisWidget->SetViewport( 0.0, 0.0, 0.15, 0.15 );
	Pool::mainAxisWidget->SetEnabled( 1 );
	Pool::mainAxisWidget->InteractiveOn();

	/******************DISPLAY DISTANCE**************************/
	
	Pool::distanceWidget->SetInteractor(Pool::interactor);
	Pool::distanceWidget->SetRepresentation(Pool::rep3D);
	Pool::distanceWidget->CreateDefaultRepresentation();
	Pool::distanceWidget->On();

	/*********INITIALIZING COVER BALLS***********/
	for (int i=0; i<Pool::coverBallActors.size(); i++) {
		//cout<<"added!"<<endl;
		Pool::renderer->AddActor(Pool::coverBallActors[i]);
		Pool::renderer->AddActor(Pool::coverCenterActors[i]);
	}

	/*********INITIALIZE TEXT ACTOR***************/
	Pool::thetaActor->GetTextProperty()->SetFontSize ( 24 );
	Pool::thetaActor->SetPosition ( 10, Pool::renderWindow->GetSize()[1]-30 );
	Pool::thetaActor->GetTextProperty()->SetColor ( 1.0,0.0,0.0 );

	Pool::phiActor->GetTextProperty()->SetFontSize ( 24 );
	Pool::phiActor->SetPosition ( 10,  Pool::renderWindow->GetSize()[1]-60 );
	Pool::phiActor->GetTextProperty()->SetColor ( 1.0,0.0,0.0 );

	Pool::renderer->AddActor(Pool::thetaActor);
	Pool::renderer->AddActor(Pool::phiActor);

// 	vtkSmartPointer<vtkTextActor> textActor = 
// 		vtkSmartPointer<vtkTextActor>::New();
// 	textActor->GetTextProperty()->SetFontSize ( 24 );
// 	textActor->SetPosition2 ( 10, 40 );
// 	Pool::renderer->AddActor2D ( textActor );
// 	textActor->SetInput ( "Hello world" );
// 	textActor->GetTextProperty()->SetColor ( 1.0,0.0,0.0 );

    
    for (int i=0; i<Pool::sceneActors.size(); i++) {
        Pool::renderer->AddActor(Pool::sceneActors[i]);
    }
    
    for (int i=0; i<3; i++) {
        Pool::renderer->AddActor(Pool::bigSpheresActor[i]);
       // Pool::axisActors[i]->SetCamera(Pool::camera);
       // Pool::renderer->AddActor(Pool::axisActors[i]);
    }
    
    double* camPos = Pool::camera->GetPosition();
    
    //cout<<"cam pos : "<<camPos[0]<<","<<camPos[1]<<","<<camPos[2]<<endl;
}

int main()
{
    /***************TYPEDEFS & CONSTANTS******************/
    const Point3i center(0,0,0);
    /***************TYPEDEFS & CONSTANTS END***************/

    /****************VTK STANDARD STUFFS******************/
    vtkSmartPointer<myInteractorStyle>                 trackBall;
    /****************VTK STANDARD STUFFS******************/
    
    const string filename = "skin_mask.mhd";
    
    //CRD::drawer(filename);
    
    
    SetupWindow();
    
    //append scenes;
	AppendScene("hepar.mhd", 0.3, 1.0, 1.0 , 1.0);
	AppendScene("vessel_mask.mhd", 1.0, 1.0, 0.0 , 0.0);
	AppendScene("tumour.mhd", 0.2, 0.0, 0.0 , 1.0,true);
	AppendScene("skin_mask.mhd", 0.2, 1.0, 1.0, 1.0);
	AppendScene("bone_mask.mhd",0.2,0.5,0.5,0);

	/*AppendScene("3.mhd", 0.3, 1.0, 1.0 , 1.0);
	AppendScene("result_modifies_mask.mhd", 1.0, 1.0, 0.0 , 0.0);
	AppendScene("3_tumor.mhd", 0.2, 0.0, 0.0 , 1.0,true);
	AppendScene("skin_mask.mhd", 0.2, 1.0, 1.0, 1.0);*/

    SetupAABB();
    
    Pool::RefreshScene();
   
    
    Pool::interactor->Initialize();
    //renderer->Render();
    Pool::renderWindow->Render();
    Pool::interactor->Start();
    
    return EXIT_SUCCESS;
}