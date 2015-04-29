#include "Cover.h"
#include "Miniball.hpp"
#include <iostream>
#include"connectivity_tumor.h"
namespace YZX{

	void balls_Execute(
		_In_ const string& tumor_mhd_file_name,//进入的肿瘤文件名称
		_In_ const double needle_radius,
		_Out_ int& min_balls_num,//输出的最小覆盖球个数
		_Out_ CenterType& centers,//输出的球中心
		_Out_ RadiusType& radiuses//输出的球半径
		)
	{

		min_balls_num = 0;
		//读取文件并且把点存到PointSeries中去
		
		vtkMetaImageReader * reader = vtkMetaImageReader::New();
		reader->SetFileName (tumor_mhd_file_name.c_str());
		reader->Update();

		vtkImageData * data = vtkImageData::New();
		data = reader->GetOutput();

		int dims[3];

		reader->GetOutput()->GetDimensions(dims);

		Mat* Mask = new Mat[dims[2]];
		for (int i=0;i<dims[2];i++)
		{
			Mask[i].create(dims[1],dims[0],CV_8UC1);
			Mask[i].setTo(Scalar(0));
		}

		//int aa=0;
		for(int k=0; k<dims[2]; k++)

		{

			for(int j=0; j<dims[1]; j++)

			{

				for(int i=0; i<dims[0]; i++)

				{
					unsigned char * voxel =(unsigned char *) (reader->GetOutput()->GetScalarPointer(i, j, k) );
					if(  *voxel > 0)
					{
						Mask[k].at<uchar>(j,i)=1;
					   // aa++;
					}
					
					
				}

			}

		}
	//std::cout<<aa<<std::endl;
		
		vector<vector<Point3i> > connectedMap;
		connectivity::ComputeConnectedParts( Mask,dims[0], dims[1], dims[2], connectedMap );

		
		vector<Mat> PointSeries(connectedMap.size());
		vector<vtkPoints *> point(connectedMap.size());
		vector<vtkCellArray *>conn(connectedMap.size());
// 		vector<vtkPolyData *> poly(connectedMap.size());
// 		vector<vtkPolyDataMapper *> mapper1(connectedMap.size());
// 		vector<vtkActor *>actor1(connectedMap.size());
// 		vector<vector<vtkSphereSource *>> sphere(connectedMap.size());
// 		vector<vector<vtkPolyDataMapper*>>mapper(connectedMap.size());
// 		vector<vector<vtkActor*>>actor(connectedMap.size());
// 		vtkRenderer * renderer = vtkRenderer::New();
		for(unsigned int i=0;i<connectedMap.size();i++)
		{
			

		PointSeries[i].create(connectedMap[i].size(),3,CV_32FC1);

		 point[i] =vtkPoints::New();
	
		int aaaa=0;

		for(unsigned int j=0;j<connectedMap[i].size();j++)
		{
			
			double p[3];
			data->GetPoint(connectedMap[i][j].z*dims[1]*dims[0]+connectedMap[i][j].y*dims[0]+connectedMap[i][j].x,p);

			for (int dim=0;dim<3;dim++)
			{
				PointSeries[i].at<float>(aaaa,dim) = p[dim];

			}
			++aaaa;
			point[i]->InsertNextPoint(p[0],p[1],p[2]);
		
		}
		bool isSatisfied=false;
		int  ball_num = 1;
	
		vector<double> radius;
		vector<vector<double> > out_centers;

		

		while(!isSatisfied)
		{
			double maxR = -1;
			
			GetBallCover(PointSeries[i],ball_num,out_centers,radius,maxR);

			//centers.clear();
			//centers.resize(ball_num);
			//radiuses.clear();
			//radiuses.resize(ball_num);
	
// 			radius.clear();
// 			out_centers.clear();
// 			radius.resize(ball_num);
// 			out_centers.resize(ball_num);

			if(maxR<needle_radius)
			{
				
				min_balls_num += ball_num;
				centers.insert(centers.end(),out_centers.begin(),out_centers.end());
				radiuses.insert(radiuses.end(),radius.begin(),radius.end());

// 				for(int ii=0;ii<min_balls_num;ii++)
// 				{
// 					centers[ii].resize(3,0);
// 					for(int j=0;j<3;j++)
// 					{
// 						
// 						centers[ii][j] = out_centers[ii][j];
// 						
// 					}
// 				
// 					radiuses[ii] = radius[ii];
// 				
//				}
				break;
				
			}

			ball_num++;

		}

		/*std::cout<<"最小覆盖球个数"<<min_balls_num<<::endl;
		for(int iii=0;iii<min_balls_num;iii++)
		{
		std::cout<<"球"<<iii+1<<"半径：   "<<radiuses[iii]<<";   球心坐标：  ("<<centers[iii][0]<<","<<centers[iii][1]<<","<<centers[iii][2]<<")"<<std::endl;

		}*/
#if 0
		conn[i] =vtkCellArray::New();
		for(unsigned int k=0;k<connectedMap[i].size();k++)
		{
			vtkIdType id = k;
			conn[i]->InsertNextCell(1,&id);
		}

		

		 poly[i] = vtkPolyData::New();
		poly[i]->SetPoints(point[i]);
		poly[i]->SetVerts(conn[i]);

		 mapper1[i] = vtkPolyDataMapper::New();
		mapper1[i]->SetInput(poly[i]);

		actor1[i] = vtkActor::New();
		actor1[i]->SetMapper(mapper1[i]);



		 sphere[i].resize(min_balls_num);   
		for(int j=0;j<min_balls_num;j++)
		{
			sphere[i][j]=vtkSphereSource::New();
			sphere[i][j]->SetThetaResolution(100);
			sphere[i][j]->SetPhiResolution(50);
			sphere[i][j]->SetRadius(radiuses[j]);
			sphere[i][j]->SetCenter(centers[j][0],centers[j][1],centers[j][2]);
		}



	     mapper[i].resize(min_balls_num);
		for(int j=0;j<min_balls_num;j++)
		{
			mapper[i][j]=vtkPolyDataMapper::New();
			mapper[i][j]->SetInputConnection(sphere[i][j]->GetOutputPort());
		}


		actor[i].resize(min_balls_num);
		for(int j=0;j<min_balls_num;j++)
		{
			actor[i][j]=vtkActor::New();
			actor[i][j]->SetMapper(mapper[i][j]);
			actor[i][j]->GetProperty()->SetOpacity(0.7);
			actor[i][j]->GetProperty()->SetColor(1,0,0);
		}

		

		
		for(int j=0;j<min_balls_num;j++)
		{
			renderer->AddActor(actor[i][j]);
			renderer->AddActor(actor1[i]);
			
		}

	}

		vtkRenderWindow * renWin = vtkRenderWindow::New();
		renWin->AddRenderer(renderer);
		renWin->SetSize(500,500);
		renWin->Render();

		vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
		iren->SetRenderWindow(renWin);
		iren->Initialize();
		iren->Start();
#endif
		}


     
	}

	void GetBallCover(
		_In_ const Mat& PointSeries,
		_In_ const int clusterNum,
		_Out_ CenterType& out_centers,
		_Out_ RadiusType& out_radiuses,
		_Out_ double& maxR
		)
	{

		out_centers.clear();
		out_centers.resize(clusterNum);
		out_radiuses.clear();
		out_radiuses.resize(clusterNum);
		const int n = PointSeries.rows;      // number of points

	

		typedef double mytype;            // coordinate type

		Mat labels;

		//cluster points according to euclidean distance
		kmeans(PointSeries, clusterNum, labels,
			TermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 100, 1e-10),
			10, KMEANS_PP_CENTERS);

		//draw these points on an image
		vector<int> labelCount(clusterNum,0);
		for (int i=0;i<n;i++)
		{
			int label = labels.at<int>(i,0);
			labelCount[label]++;
		}

		// define the types of iterators through the points and their coordinates
		// ----------------------------------------------------------------------
		typedef mytype* const* PointIterator; 
		typedef const mytype* CoordIterator;

		typedef Miniball::Miniball<Miniball::CoordAccessor<PointIterator, CoordIterator> > MB;

		vector<double*> centers(clusterNum,NULL);
		vector<double**> cPoints(clusterNum,NULL);

		for (int i=0;i<clusterNum;i++)
		{
			cPoints[i]=new double*[labelCount[i]];
			for (int j=0;j<labelCount[i];j++)
			{
				cPoints[i][j] = new double[3];
			}
		}

		vector<int> counters(clusterNum,0);

		for (int i=0;i<n;i++)
		{
			int label = labels.at<int>(i,0);

			for (int dim=0;dim<3;dim++)
			{
				cPoints[label][counters[label]][dim]=PointSeries.at<float>(i,dim);
			}

			counters[label]++;
		}

		double maxRadius = -1;
		for (int i=0;i<clusterNum;i++)
		{

			MB mb (3, cPoints[i], cPoints[i]+labelCount[i]);
			double radius = sqrt(mb.squared_radius());
			if(radius>maxRadius) maxRadius=radius;

			out_radiuses[i] = radius;
			const double* current_center = mb.center();
			out_centers[i].resize(3,0);
			for (int j=0;j<3;j++)
			{
				out_centers[i][j] = current_center[j];
			}
		}

	

		maxR = maxRadius;


	}

}