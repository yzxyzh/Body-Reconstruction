#include"drawing.h"

#include <vtkSmartPointer.h>
#include "Pool.h"

namespace CRD{

	//void drawer(_In_ const vector<string>& fileNames)
	void drawer(_In_ const string& skin_mhd_file_name)
	{

		vector<vtkSmartPointer<vtkActor> > polyActors;
		vector<vtkSmartPointer<vtkPolyData> > polyDatas;




		vtkMetaImageReader* reader =

				vtkMetaImageReader::New();

		reader->SetFileName(skin_mhd_file_name.c_str());

		reader->Update();


		int dim[3];

		int extent[6];


		double spacing[3];


		double origin[3];


		reader->GetOutput()->GetDimensions(dim);
		//cout<<"intput dimension"<<dim[0]<<" "<<dim[1]<<" "<<dim[2]<<endl;

		reader->GetOutput()->GetExtent(extent);
		//cout<<"input extension"<<extent[0]<<" "<<extent[1]<<" "<<extent[2]<<" "<<extent[3]<<" "<<extent[4]<<" "<<extent[5]<<endl;

		reader->GetOutput()->GetSpacing(spacing);
		//cout<<"input spacing"<<spacing[0]<<" "<<spacing[1]<<" "<<spacing[2]<<endl;

		reader->GetOutput()->GetOrigin(origin);
		//cout<<"input origin"<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<endl;


		double center[3];

		center[0] = origin[0] + spacing[0] * 0.5 * (extent[0] + extent[1]);

		center[1] = origin[1] + spacing[1] * 0.5 * (extent[2] + extent[3]);

		center[2] = origin[2] + spacing[2] * 0.5 * (extent[4] + extent[5]);

		//cout<<"sphere center"<<center[0]<<" "<<center[1]<<" "<<center[2]<<endl;


		double xmin=-1,xmax=1,ymin=-1,ymax=1,zmin=-1,zmax=1;

		vtkImageData * image = vtkImageData::New();
		image=reader->GetOutput();


		vtkMarchingCubes * iso = vtkMarchingCubes::New();
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


		vtkSmoothPolyDataFilter* smoothFilter =
				vtkSmoothPolyDataFilter::New();
		smoothFilter->SetInputConnection(iso->GetOutputPort());
		smoothFilter->SetNumberOfIterations(10);
		smoothFilter->FeatureEdgeSmoothingOn();
		smoothFilter->BoundarySmoothingOn();
		smoothFilter->Update();

		vtkQuadricClustering *decimate =vtkQuadricClustering::New();
		decimate->SetNumberOfXDivisions(1200);
		decimate->SetNumberOfYDivisions(1200);
		decimate->SetNumberOfZDivisions(1200);

#if VTK_MAJOR_VERSION <6
    decimate->SetInput(smoothFilter->GetOutput());
#else
		decimate->SetInputData(smoothFilter->GetOutput());
#endif



		vtkPolyData * data = vtkPolyData::New();//vtkPolyDataµƒΩ”ø⁄°™°™°™°™CRD
		data=smoothFilter->GetOutput();

		//double xmin=0,xmax=0,ymax=0,ymin=0,zmax=0,zmin=0;

		double* bounds = data->GetBounds();

		/*std::cout  << "xmin: " << bounds[0] << " "
		<< "xmax: " << bounds[1] << std::endl
		<< "ymin: " << bounds[2] << " "
		<< "ymax: " << bounds[3] << std::endl
		<< "zmin: " << bounds[4] << " "
		<< "zmax: " << bounds[5] << std::endl;
*/
		xmin = bounds[0];
		xmax = bounds[1];
		ymin = bounds[2];
		ymax = bounds[3];
		zmin = bounds[4];
		zmax = bounds[5];
		double xcenter,ycenter,zcenter,radius;
		xcenter=0.5*(xmax+xmin);
		ycenter=0.5*(ymax+ymin);
		zcenter=0.5*(zmax+zmin);
		radius=0.5*sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin));



		vtkRegularPolygonSource * polygonSource1 =vtkRegularPolygonSource::New();//this is used for circle drawing--CRD

		polygonSource1->GeneratePolygonOff();
		polygonSource1->SetNumberOfSides(50);
		polygonSource1->SetRadius(radius);
		polygonSource1->SetCenter(xcenter,ycenter,zcenter);
		polygonSource1->SetNormal(0,0,1);
		polygonSource1->Update();

		vtkRegularPolygonSource * polygonSource2 =vtkRegularPolygonSource::New();//this is used for circle drawing--CRD

		polygonSource2->GeneratePolygonOff();
		polygonSource2->SetNumberOfSides(50);
		polygonSource2->SetRadius(radius);
		polygonSource2->SetCenter(xcenter,ycenter,zcenter);
		polygonSource2->SetNormal(0,1,0);
		polygonSource2->Update();

		vtkRegularPolygonSource * polygonSource3 =vtkRegularPolygonSource::New();//this is used for circle drawing--CRD

		polygonSource3->GeneratePolygonOff();
		polygonSource3->SetNumberOfSides(50);
		polygonSource3->SetRadius(radius);
		polygonSource3->SetCenter(xcenter,ycenter,zcenter);
		polygonSource3->SetNormal(1,0,0);
		polygonSource3->Update();


		/*  vtkVectorNorm * gradient = vtkVectorNorm::New();
           gradient->SetInputConnection(decimate->GetOutputPort());*/


		/*vtkPolyDataNormals * normal = vtkPolyDataNormals::New();
           normal->SetInputConnection(smoothFilter->GetOutputPort());
           normal->SetFeatureAngle(60.0);*/

		/*vtkDataSetMapper * cubeMapper = vtkDataSetMapper::New();
           cubeMapper->SetInputConnection(decimate->GetOutputPort());*/

		vtkPolyDataMapper * dataMapper = vtkPolyDataMapper::New();
		dataMapper->SetInputConnection(decimate->GetOutputPort());

		vtkPolyDataMapper * dataMapper1 = vtkPolyDataMapper::New();
		dataMapper1->SetInputConnection(polygonSource1->GetOutputPort());

		vtkPolyDataMapper * dataMapper2 = vtkPolyDataMapper::New();
		dataMapper2->SetInputConnection(polygonSource2->GetOutputPort());

		vtkPolyDataMapper * dataMapper3 = vtkPolyDataMapper::New();
		dataMapper3->SetInputConnection(polygonSource3->GetOutputPort());



		vtkActor * cubeActor = vtkActor::New();
		cubeActor->SetMapper(dataMapper);
		cubeActor->GetProperty()->SetOpacity(0.9);
		cubeActor->GetProperty()->SetColor(0.5,0.5,0.5);

		vtkActor * cubeActor1 = vtkActor::New();
		cubeActor1->SetMapper(dataMapper1);

		vtkActor * cubeActor2 = vtkActor::New();
		cubeActor2->SetMapper(dataMapper2);

		vtkActor * cubeActor3 = vtkActor::New();
		cubeActor3->SetMapper(dataMapper3);

		vtkRenderer* renderer =

				vtkRenderer::New();

		renderer->AddActor(cubeActor);
		renderer->AddActor(cubeActor1);
		renderer->AddActor(cubeActor2);
		renderer->AddActor(cubeActor3);
		renderer->SetBackground(0.0,0.0,0.0);



		vtkRenderWindow * renderWindow =

				vtkRenderWindow::New();

		renderWindow->SetSize(600, 600);

		renderWindow->AddRenderer(renderer);

		renderWindow->Render();


		vtkRenderWindowInteractor * renderWindowInteractor =

				vtkRenderWindowInteractor::New();


		renderWindowInteractor->SetRenderWindow(renderWindow);



		vtkCamera * aCamera = vtkCamera::New();
		aCamera->ComputeViewPlaneNormal();
		aCamera->ParallelProjectionOn();

		renderer->SetActiveCamera(aCamera);
		renderer->ResetCamera();


		renderWindowInteractor->Initialize();
		renderWindowInteractor->Start();


	}
}