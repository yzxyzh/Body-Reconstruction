#include "Pool.h"

#include <vtkProperty.h>
#include <string>
#include <sstream>

using namespace std;

std::string any2Str(double in)//any2Str function  Converts an anytype value to a str value
{
	stringstream ss;
	ss<<in;
	std::string str;
	ss>>str;
	return str;
	
}


vtkSmartPointer<vtkRenderer>                       Pool::renderer =
vtkSmartPointer<vtkRenderer>::New();
vtkSmartPointer<vtkRenderWindow>                   Pool::renderWindow =
vtkSmartPointer<vtkRenderWindow>::New();
vtkSmartPointer<vtkRenderWindowInteractor>         Pool::interactor =
vtkSmartPointer<vtkRenderWindowInteractor>::New();
vtkSmartPointer<vtkCamera>                         Pool::camera =
vtkSmartPointer<vtkCamera>::New();
vtkSmartPointer<vtkSphereSource>                   Pool::sphere =
vtkSmartPointer<vtkSphereSource>::New();
vtkSmartPointer<vtkActor>                          Pool::sphereActor =
vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkPolyDataMapper>                 Pool::sphereMapper =
vtkSmartPointer<vtkPolyDataMapper>::New();
vtkSmartPointer<vtkSphereSource>                   Pool::Point =
vtkSmartPointer<vtkSphereSource>::New();
vtkSmartPointer<vtkActor>                          Pool::PointActor =
vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkPolyDataMapper>                 Pool::PointMapper =
vtkSmartPointer<vtkPolyDataMapper>::New();
vtkSmartPointer<vtkTransform>                      Pool::userTransform =
vtkSmartPointer<vtkTransform>::New();
vtkSmartPointer<vtkCylinderSource>                 Pool::cylinder =
vtkSmartPointer<vtkCylinderSource>::New();
vtkSmartPointer<vtkPolyDataMapper>                 Pool::cylinderMapper =
vtkSmartPointer<vtkPolyDataMapper>::New();
vtkSmartPointer<vtkActor>                          Pool::cylinderActor =
vtkSmartPointer<vtkActor>::New();
vtkMatrix4x4*                                      Pool::cameraTransform =
vtkMatrix4x4::New();
vtkSmartPointer<vtkLineSource>                     Pool::line =
vtkSmartPointer<vtkLineSource>::New();
vtkSmartPointer<vtkPolyDataMapper>                 Pool::lineMapper =
vtkSmartPointer<vtkPolyDataMapper>::New();
vtkSmartPointer<vtkActor>                          Pool::lineActor =
vtkSmartPointer<vtkActor>::New();
vector<vtkSmartPointer<vtkActor> >                 Pool::sceneActors =
vector<vtkSmartPointer<vtkActor> >();
vector<vtkSmartPointer<vtkPolyData> >              Pool::scenePolyDatas =
vector<vtkSmartPointer<vtkPolyData> > ();
vector<vtkSmartPointer<vtkRegularPolygonSource> >  Pool::bigSpheres =
vector<vtkSmartPointer<vtkRegularPolygonSource> >(3);
vector<vtkSmartPointer<vtkPolyDataMapper> >        Pool::bigSpheresMapper =
vector<vtkSmartPointer<vtkPolyDataMapper> >(3);
vector<vtkSmartPointer<vtkActor> >                 Pool::bigSpheresActor =
vector<vtkSmartPointer<vtkActor> >(3);
vector<vtkSmartPointer<vtkAxisActor> >             Pool::axisActors =
vector<vtkSmartPointer<vtkAxisActor> >(3);
vtkSmartPointer<vtkAnnotatedCubeActor>             Pool::mainAxisActor =
vtkSmartPointer<vtkAnnotatedCubeActor>::New();
vtkSmartPointer<vtkOrientationMarkerWidget>        Pool::mainAxisWidget = 
vtkSmartPointer<vtkOrientationMarkerWidget>::New();
vtkSmartPointer<vtkDistanceWidget>                 Pool::distanceWidget=
vtkSmartPointer<vtkDistanceWidget>::New();
vtkSmartPointer<vtkDistanceRepresentation3D>       Pool::rep3D=
vtkSmartPointer<vtkDistanceRepresentation3D>::New();
//cover balls;
vector<vtkSmartPointer<vtkSphereSource> >          Pool::coverBallSources =
	vector<vtkSmartPointer<vtkSphereSource> >();
vector<vtkSmartPointer<vtkPolyDataMapper> >        Pool::coverBallMappers =
	vector<vtkSmartPointer<vtkPolyDataMapper> >();
vector<vtkSmartPointer<vtkActor> >                 Pool::coverBallActors =
	vector<vtkSmartPointer<vtkActor> >();
vtkSmartPointer<vtkTextActor>                      Pool::thetaActor = 
	vtkSmartPointer<vtkTextActor>::New();
vtkSmartPointer<vtkTextActor>                      Pool::phiActor = 
vtkSmartPointer<vtkTextActor>::New();

vector<vtkSmartPointer<vtkSphereSource> >          Pool::coverCenterSources =
	vector<vtkSmartPointer<vtkSphereSource> >();
vector<vtkSmartPointer<vtkPolyDataMapper> >        Pool::coverCenterMappers =
	vector<vtkSmartPointer<vtkPolyDataMapper> >();
vector<vtkSmartPointer<vtkActor> >                 Pool::coverCenterActors =
	vector<vtkSmartPointer<vtkActor> >();

int                                                Pool::currentBallIndex = 
0;

vector<Point3d> Pool::tmpCenter = vector<Point3d>();
int Pool::ballNumber=0;


myInteractorStyle* Pool::trackball = myInteractorStyle::New();

Point3d Pool::currPoint = Point3d(0,0,1);
Point3d Pool::sphereCenter = Point3d(0,0,0);

Point3d Pool::tumorCenter = Point3d(0,0.3,0);

vector<vtkSmartPointer<vtkActor> > Pool::lastIntersectionActors = vector<vtkSmartPointer<vtkActor> >();

vector<vector<Point3d> > Pool::interSectors = vector<vector<Point3d> >();
vector<vector<double> > Pool::colors = vector<vector<double> >();

void Pool::RefreshScene()
{
    renderer->RemoveActor(Pool::PointActor);
    for (int i=0; i<lastIntersectionActors.size(); i++) {
        renderer->RemoveActor(lastIntersectionActors[i]);
    }
    renderer->RemoveActor(Pool::lineActor);

	renderer->RemoveActor(thetaActor);
	renderer->RemoveActor(phiActor);

    lastIntersectionActors.clear();
    Pool::Point->SetCenter(Pool::currPoint.x, Pool::currPoint.y, Pool::currPoint.z);
    Pool::Point->SetRadius(2.0);
    Pool::PointMapper->SetInputConnection(Pool::Point->GetOutputPort());
    Pool::PointActor->SetMapper(Pool::PointMapper);
    Pool::PointActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
    
    
    Pool::line->SetPoint1(Pool::currPoint.x, Pool::currPoint.y, Pool::currPoint.z);
    Pool::line->SetPoint2(Pool::tumorCenter.x, Pool::tumorCenter.y, Pool::tumorCenter.z);
    Pool::lineMapper->SetInputConnection(Pool::line->GetOutputPort());
    Pool::lineActor->SetMapper(Pool::lineMapper);
    Pool::lineActor->GetProperty()->SetColor(0.0, 0.0, 1.0);
    Pool::lineActor->GetProperty()->SetLineWidth(3.0);
    
    renderer->AddActor(lineActor);
    //update intersector points
    int size = Pool::interSectors.size();
    
    for (int i=0; i<size; i++) {
        
        for (int j=0; j<Pool::interSectors[i].size(); j++) {
            
            Point3d currPoint = Pool::interSectors[i][j];
            vtkSmartPointer<vtkSphereSource>   sphere = vtkSmartPointer<vtkSphereSource>::New();
            vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            vtkSmartPointer<vtkActor>          sphereActor  = vtkSmartPointer<vtkActor>::New();
            
            sphere->SetCenter(currPoint.x, currPoint.y, currPoint.z);
            sphere->SetRadius(3.0);
            
            sphereMapper->SetInputConnection(sphere->GetOutputPort());
            
            sphereActor->SetMapper(sphereMapper);
            sphereActor->GetProperty()->SetColor(Pool::colors[i][0], Pool::colors[i][1], Pool::colors[i][2]);
            
            lastIntersectionActors.push_back(sphereActor);
            
            renderer->AddActor(sphereActor);
    
        }
        
    }
    
    
    renderer->AddActor(PointActor);

	//update theta & phi
	//theta:angle between line and x-z plane
	//phi:angle between projected line and positiv x axis;
	string thetaStr = "THETA angle = ";
	string phiStr = "PHI angle = ";
	
	double theta,phi;

	Point3d vec = currPoint-tumorCenter;

	//calculate theta first, theta has range [-PI/2,PI/2];
	double r = cv::norm(vec);
	double arcsin = vec.y/r;
	double angle = asin(arcsin);
	theta = angle*180/CV_PI;

	//then calcualte phi,range from [0,2PI];
	double r_projection = sqrt(vec.x*vec.x + vec.z*vec.z);
	if(vec.z>=0)
	{ 
		double arccos = vec.x/r_projection;
		double angle = acos(arccos);
		phi = angle*180/CV_PI;
	}else{
		double arccos = -vec.x/r_projection;
		double angle = acos(arccos);
		angle+=CV_PI;
		phi = angle*180/CV_PI;
	}

	thetaStr+=any2Str(theta);
	phiStr+=any2Str(phi);

	thetaActor->SetInput(thetaStr.c_str());
	phiActor->SetInput(phiStr.c_str());

	renderer->AddActor(thetaActor);
	renderer->AddActor(phiActor);


}









