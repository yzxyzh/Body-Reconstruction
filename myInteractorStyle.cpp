#include "myInteractorStyle.h"

#include "Pool.h"
#include <vtkMatrix4x4.h>

vtkStandardNewMacro(myInteractorStyle);

myInteractorStyle::myInteractorStyle():
vtkInteractorStyleTrackballCamera()
{
    this->isMovingPoint = false;
    this->currentPoint = Point3d(0,0,1);
    this->radius = 1.0;
    this->pointMoving = false;
    this->theta = 0.0;
    this->phi = 0.0;
    sphereActor = NULL;
    this->tmpCamera = vtkCamera::New();
}

void myInteractorStyle::UpdateIntersectors()
{
    int obj_size = octTreeList.size();
    
    Pool::interSectors.clear();
    
    if(obj_size == 0) return;
    
    Pool::interSectors.resize(obj_size);
    
    double curr[3]={Pool::currPoint.x,Pool::currPoint.y,Pool::currPoint.z};
    double dst[3]={Pool::tumorCenter.x,Pool::tumorCenter.y,Pool::tumorCenter.z};
    
    vtkSmartPointer<vtkPoints> intPoint = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    
    for (int k=0; k<obj_size; k++) {
        
        vtkSmartPointer<vtkOBBTree> tree = octTreeList[k];
        
        tree->IntersectWithLine(curr, dst, intPoint, cellIds);
        
        int pSize = intPoint->GetNumberOfPoints();
        
        for (int i=0; i<pSize; i++) {
            
            vtkIdType id = i;
            
            double* currIntP = intPoint->GetPoint(id);
            
            Point3d newP(currIntP[0],currIntP[1],currIntP[2]);
            
            Pool::interSectors[k].push_back(newP);
            
        }
        
    }
    
   

    
}

void myInteractorStyle::AddObject(vtkPolyData *polyData)
{
    vtkSmartPointer<vtkOBBTree> tree = vtkSmartPointer<vtkOBBTree>::New();
    tree->SetDataSet(polyData);
    tree->BuildLocator();
    
    octTreeList.push_back(tree);
}

myInteractorStyle::~myInteractorStyle()
{
    if(this->tmpCamera) this->tmpCamera->Delete();
}

void myInteractorStyle::OnRightButtonDown()
{
    isMovingPoint = !isMovingPoint;
    pointMoving = false;
}



void myInteractorStyle::GetCameraTransform()
{
    vtkMatrix4x4* mat = Pool::camera->GetViewTransformMatrix();
    Pool::cameraTransform->DeepCopy(mat);
    Pool::cameraTransform->SetElement(0, 3, 0.0);
    Pool::cameraTransform->SetElement(1, 3, 0.0);
    Pool::cameraTransform->SetElement(2, 3, 0.0);

}

//void myInteractorStyle::OnMiddleButtonDown()
//{
//    this->FindPokedRenderer(this->Interactor->GetEventPosition()[0],
//                            this->Interactor->GetEventPosition()[1]);
//    if (this->CurrentRenderer == NULL)
//    {
//        return;
//    }
//    
//    this->GrabFocus(this->EventCallbackCommand);
//    this->StartPan();
//}

void myInteractorStyle::OnMouseMove()
{
    //cout<<this->State<<" PAN "<<VTKIS_PAN<<endl;
    if(!isMovingPoint || !pointMoving)
    {
        //cout<<"mouse move!"<<endl;
        
        superclass::OnMouseMove();

        //cout<<"actual Position : "<<Pool::currPoint.x<<" "<<Pool::currPoint.y<<" "<<Pool::currPoint.z<<endl;
        
    }else{
        
        //cout<<"basePoint = "<<basePoint[0]<<" "<<basePoint[1]<<" "<<basePoint[2]<<endl;
        
        vtkRenderWindowInteractor *rwi = this->Interactor;
        
        int dx = rwi->GetEventPosition()[0] - rwi->GetLastEventPosition()[0];
        int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];
        
        //cout<<"dx = "<<dx<<" dy = "<<dy<<endl;
        
        int *size = Pool::renderWindow->GetSize();
        
        //cout<<"size = "<<size[0]<<"x"<<size[1]<<endl;
        
        double delta_elevation = -20.0 / size[1];
        double delta_azimuth = -20.0 / size[0];
        
        double rxf = dx * delta_azimuth * this->MotionFactor;
        double ryf = dy * delta_elevation * this->MotionFactor;
        
        
        //tmpCamera->DeepCopy(Pool::camera);
       
        tmpCamera->Azimuth(rxf);
        tmpCamera->Elevation(ryf);
        tmpCamera->OrthogonalizeViewUp();
        
        //Pool::camera->Print(cout);
        
//        double *center = this->CurrentRenderer->GetCenter();
//        
//        double newAngle =
//        vtkMath::DegreesFromRadians( atan2( rwi->GetEventPosition()[1] - center[1],
//                                           rwi->GetEventPosition()[0] - center[0] ) );
//        
//        double oldAngle =
//        vtkMath::DegreesFromRadians( atan2( rwi->GetLastEventPosition()[1] - center[1],
//                                           rwi->GetLastEventPosition()[0] - center[0] ) );
//        
//        //vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
//        tmpCamera->Roll( newAngle - oldAngle );
//        tmpCamera->OrthogonalizeViewUp();

        
        vtkMatrix4x4* mat1 = tmpCamera->GetViewTransformMatrix();
        vtkMatrix4x4* mat2 = vtkMatrix4x4::New();
        mat2->DeepCopy(mat1);
        mat2->SetElement(0, 3, 0.0);
        mat2->SetElement(1, 3, 0.0);
        mat2->SetElement(2, 3, 0.0);
     
        //mat2->Print(cout);
        //Pool::cameraTransform->Print(cout);
        
        
        double camera_pos[4];
        double actual_pos[4];
        mat2->MultiplyPoint(basePoint, camera_pos);
        
        vtkMatrix4x4* camInv = vtkMatrix4x4::New();
        vtkMatrix4x4::Invert(Pool::cameraTransform, camInv);
        camInv->MultiplyPoint(camera_pos, actual_pos);
        
        
        //Pool::cameraTransform->MultiplyPoint(camera_pos, actual_pos);
        
        
        
        Pool::currPoint.x=actual_pos[0]+Pool::sphereCenter.x;
        Pool::currPoint.y=actual_pos[1]+Pool::sphereCenter.y;
        Pool::currPoint.z=actual_pos[2]+Pool::sphereCenter.z;
        
        UpdateIntersectors();
        
        Pool::RefreshScene();
        mat2->Delete();
        camInv->Delete();
        rwi->Render();
        
        //cout<<"did this!"<<endl;
    }
    
}

 void myInteractorStyle::OnChar()
{
     switch(Pool::interactor->GetKeyCode())
	 {
	 case 't':
	 case 'T':
		{
			Pool::currentBallIndex = (1+Pool::currentBallIndex)% (Pool::tmpCenter.size());


		 Pool::tumorCenter.x = Pool::tmpCenter[Pool::currentBallIndex].x;
		 Pool::tumorCenter.y = Pool::tmpCenter[Pool::currentBallIndex].y;
		 Pool::tumorCenter.z = Pool::tmpCenter[Pool::currentBallIndex].z;
		 
		 UpdateIntersectors();

		 Pool::RefreshScene();
		 Pool::renderWindow->Render();
		 
		}
		 break;
	 default:
		 superclass::OnChar();
		 break;
	 }

}


void myInteractorStyle::OnLeftButtonUp()
{
    if(!isMovingPoint || !pointMoving)
    {
        superclass::OnLeftButtonUp();
        
//        double p[4];
//        p[0] = Pool::currPoint.x-Pool::sphereCenter.x;
//        p[1] = Pool::currPoint.y-Pool::sphereCenter.y;
//        p[2] = Pool::currPoint.z-Pool::sphereCenter.z;
//        p[3] = 1;
//        
//        
//        double actual_pos[4];
//        GetCameraTransform();
//        Pool::cameraTransform->MultiplyPoint(p, actual_pos);
//        
//        actual_pos[0]+=Pool::sphereCenter.x;
//        actual_pos[1]+=Pool::sphereCenter.y;
//        actual_pos[2]+=Pool::sphereCenter.z;
//        
//        Pool::currPoint.x = actual_pos[0];
//        Pool::currPoint.y = actual_pos[1];
//        Pool::currPoint.z = actual_pos[2];
//
//        Pool::camera->SetFocalPoint(Pool::sphereCenter.x, Pool::sphereCenter.y, Pool::sphereCenter.z);
//        Pool::camera->SetViewUp(0, 1, 0);
//        Pool::camera->SetPosition(0, 0, 10);
//        Pool::RefreshScene();
        
        
    }else{
        //Pool::camera->DeepCopy(tmpCamera);
        pointMoving = false;
    }
    
    
}

void myInteractorStyle::OnLeftButtonDown()
{
    
    if(!isMovingPoint)
    {
        superclass::OnLeftButtonDown();
        

    }else{
        //cout<<"lB Down!"<<endl;
//        
//       
//        double p[4];
//        p[0] = Pool::currPoint.x-Pool::sphereCenter.x;
//        p[1] = Pool::currPoint.y-Pool::sphereCenter.y;
//        p[2] = Pool::currPoint.z-Pool::sphereCenter.z;
//        p[3] = 1;
//        
//        
//        double actual_pos[4];
//        GetCameraTransform();
//        Pool::cameraTransform->MultiplyPoint(p, actual_pos);
//        
//        actual_pos[0]+=Pool::sphereCenter.x;
//        actual_pos[1]+=Pool::sphereCenter.y;
//        actual_pos[2]+=Pool::sphereCenter.z;
//        
//        Pool::currPoint.x = actual_pos[0];
//        Pool::currPoint.y = actual_pos[1];
//        Pool::currPoint.z = actual_pos[2];
//        
//
//        Pool::renderer->ResetCamera();

//        double p[4];
//        p[0] = Pool::currPoint.x-Pool::sphereCenter.x;
//        p[1] = Pool::currPoint.y-Pool::sphereCenter.y;
//        p[2] = Pool::currPoint.z-Pool::sphereCenter.z;
//        p[3] = 1;
//        
//        
//        double actual_pos[4];
//        GetCameraTransform();
//        Pool::cameraTransform->MultiplyPoint(p, actual_pos);
//        
//        actual_pos[0]+=Pool::sphereCenter.x;
//        actual_pos[1]+=Pool::sphereCenter.y;
//        actual_pos[2]+=Pool::sphereCenter.z;
        
        //Pool::currPoint.x = actual_pos[0];
        //Pool::currPoint.y = actual_pos[1];
        //Pool::currPoint.z = actual_pos[2];

        
        basePoint[0] = Pool::currPoint.x-Pool::sphereCenter.x;
        basePoint[1] = Pool::currPoint.y-Pool::sphereCenter.y;
        basePoint[2] = Pool::currPoint.z-Pool::sphereCenter.z;
        basePoint[3] = 1;
        
        tmpCamera->SetPosition(Pool::camera->GetPosition());
        tmpCamera->SetViewUp(Pool::camera->GetViewUp());
        tmpCamera->SetFocalPoint(Pool::camera->GetFocalPoint());
        
        GetCameraTransform();
//        tmpCamera->DeepCopy(Pool::camera);
//        tmpCamera->Delete();
//        tmpCamera = vtkCamera::New();
//        tmpCamera->SetViewUp (0 , 1 , 0);
//        tmpCamera->SetPosition (0, 0, 10.0);
//        tmpCamera->SetFocalPoint (Pool::sphereCenter.x, Pool::sphereCenter.y, Pool::sphereCenter.z);
//        tmpCamera->ParallelProjectionOn();
        //GetCameraTransform();
        //tmpCamera->DeepCopy(Pool::camera);
        pointMoving = true;
    }
 
}