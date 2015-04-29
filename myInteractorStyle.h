#pragma once

#include <opencv2/core/core.hpp>
//#include <vtkInteractionStyleModule.h> // For export macro
#include <vtkInteractorStyle.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCamera.h>
#include <vtkCallbackCommand.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

#include <vtkActor.h>
#include <vtkSphereSource.h>
#include <vtkOBBTree.h>
#include <vtkSmartPointer.h>


using namespace cv;

class myInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    typedef vtkInteractorStyleTrackballCamera superclass;
    typedef myInteractorStyle                 self;
    static myInteractorStyle *New();
    vtkTypeMacro(myInteractorStyle,vtkInteractorStyleTrackballCamera);
    
    vtkSetMacro(isMovingPoint, bool);
    vtkSetMacro(radius, double);
    vtkSetMacro(sphereActor, vtkActor*);
    
    void OnMouseMove();
    void OnLeftButtonDown();
    void OnLeftButtonUp();
    void OnRightButtonDown();
	void OnChar();

    //void OnMiddleButtonDown();
    
    void GetCameraTransform();

    void AddObject(vtkPolyData* polyData);
    
    void UpdateIntersectors();
    
protected:
    myInteractorStyle();
    ~myInteractorStyle();
    
    Point3d currentPoint;
    bool    isMovingPoint;
    double  radius;
    
    double  theta;
    double  phi;
    
    vtkCamera*      tmpCamera;
    double   basePoint[4];
    
    bool    pointMoving;
    vtkActor* sphereActor;
    
    vector<vtkSmartPointer<vtkOBBTree> > octTreeList;
    
    
private:
    //禁止拷贝构造函数
    myInteractorStyle(const myInteractorStyle&);  // Not implemented.
    void operator=(const myInteractorStyle&);  // Not implemented.
};