#include "connectivity_tumor.h"
#include <iostream>


	namespace connectivity{

		Mat* Mask=NULL;
		int dim3=-1;
		
		int dim2=-1;
		int dim1=-1;

		typedef struct islandPoints {
			Point3i featurePoint;//�����������
			int islandVolume;//��������
			std::vector<Point3i> pointsInIsland;
		} islandPoints;

		//�õ�������ͨ��֧
		void GetIsland(Point3i& startPoint,islandPoints& island,Mat * ImageFloodFill)
		{
			island.pointsInIsland.clear();

			island.featurePoint=startPoint;
			Point3i Seed;
			Seed=startPoint;
		
			if (Mask[Seed.z].at<uchar>(Seed.y,Seed.x)==0
				|| Seed.z<0
				||Seed.z>dim3){

					std::cout<<"��ѡ��㣡"<<std::endl;
			}
			
			std::vector<Point3i> FloodFillQueue;
			Point3i pointTmp;//���ڴ������м����
			Point3i pointToPush;
			int xMax=dim1;
			int yMax=dim2;
			int zMin=0;
			int zMax=dim3;
			
			FloodFillQueue.push_back(Seed);
		

			while(!FloodFillQueue.empty())
			{
				pointTmp=FloodFillQueue.back();//�������һ�����ݣ��������������Ƿ���ڡ� 
				FloodFillQueue.pop_back();//ɾ�����һ������
				//
				if (Mask[pointTmp.z].at<uchar>(pointTmp.y,pointTmp.x)==1
					&& ImageFloodFill[pointTmp.z].at<uchar>(pointTmp.y,pointTmp.x) == 0
					){
					ImageFloodFill[pointTmp.z].at<uchar>(pointTmp.y,pointTmp.x)=1;
					island.pointsInIsland.push_back(pointTmp);
					pointToPush=pointTmp;
					pointToPush.x--;//-1
					if (pointToPush.x>=0 && 
						ImageFloodFill[pointToPush.z].at<uchar>(pointToPush.y,pointToPush.x)==0 ){
							FloodFillQueue.push_back(pointToPush);}
					pointToPush.x+=2;//+1
					if (pointToPush.x<xMax && 
						ImageFloodFill[pointToPush.z].at<uchar>(pointToPush.y,pointToPush.x)==0 ){
							FloodFillQueue.push_back(pointToPush);}
					pointToPush.x--;//��0

					pointToPush.y--;//-1
					if (pointToPush.y>=0 && 
						ImageFloodFill[pointToPush.z].at<uchar>(pointToPush.y,pointToPush.x)==0 ){
							FloodFillQueue.push_back(pointToPush);}
					pointToPush.y+=2;//+1
					if (pointToPush.y<yMax && 
						ImageFloodFill[pointToPush.z].at<uchar>(pointToPush.y,pointToPush.x)==0 ){
							FloodFillQueue.push_back(pointToPush);}
					pointToPush.y--;//��0

					pointToPush.z--;//-1
					if (pointToPush.z>=zMin && 
						ImageFloodFill[pointToPush.z].at<uchar>(pointToPush.y,pointToPush.x)==0 ){
							FloodFillQueue.push_back(pointToPush);}
					pointToPush.z+=2;//+1
					if (pointToPush.z<zMax && 
						ImageFloodFill[pointToPush.z].at<uchar>(pointToPush.y,pointToPush.x)==0 ){
							FloodFillQueue.push_back(pointToPush);}
					pointToPush.z--;//��0
					
				}
			}//end while
			island.islandVolume=island.pointsInIsland.size();
			
		}
		//�õ����е���ͨ��֧
		void ComputeConnectedParts(_in_ Mat * Mask1,_in_ int dims1,_in_ int dims2,_in_ int dims3,_out_ vector<vector<Point3i> >& connectedMap )
		{
			connectedMap.clear();
			
			//std::cout<<"hahhahah"<<dims1<<std::endl;
			dim1=dims1;
			dim2=dims2;
			dim3=dims3;//����ط�����û�б�Ҫ
			Mask=new Mat[dim3];
			for (int i=0;i<dim3;i++)
			{
				Mask[i].create(dim2,dim1,CV_8UC1);
				Mask[i].setTo(Scalar(0));
			}
			
			for(int k=0;k<dim3;k++)
			{
				for(int j=0;j<dim2;j++)
				{
					for(int i=0;i<dim1;i++)
					{
						Mask[k].at<uchar>(j,i)=Mask1[k].at<uchar>(j,i);
						
						
					}
					
				}
			}
	

			Mat* FloodFillTmp=new Mat[dim3];
			for (int i=0;i<dim3;i++)
			{
				FloodFillTmp[i].create(dim2,dim1,CV_8UC1);
				FloodFillTmp[i].setTo(Scalar(0));
			}
			

			for (int z=0;z<dim3;z++)
			{
				for (int y=0;y<dim2;y++)
				{
					for (int x=0;x<dim1;x++)
					{
						if (Mask[z].at<uchar>(y,x)>0)
						{
							
							islandPoints myIslandPoints;
							Point3i point_tmp(x,y,z);
							
							connectivity::GetIsland(point_tmp,myIslandPoints,FloodFillTmp);//ע�����(x,y)�Ǿ������꣬���ɵ���(y,x)
						
							for (int k=0;k<myIslandPoints.islandVolume;k++)
							{
								Mask[myIslandPoints.pointsInIsland[k].z].at<uchar>
									(myIslandPoints.pointsInIsland[k].y,myIslandPoints.pointsInIsland[k].x)=0;
							}
							vector<Point3i> newPart;
							newPart.clear();
							for ( unsigned int i=0;i<myIslandPoints.pointsInIsland.size();i++)
							{
								newPart.push_back(myIslandPoints.pointsInIsland[i]);
							}
							newPart.resize(newPart.size());
							connectedMap.push_back(newPart);
							for (int k=0;k<myIslandPoints.islandVolume;k++)
							{
								FloodFillTmp[myIslandPoints.pointsInIsland[k].z].at<uchar>
									(myIslandPoints.pointsInIsland[k].y,myIslandPoints.pointsInIsland[k].x)=0;
							}


						}
					}
				}
			}
			
			delete [] FloodFillTmp;
			
			std::cout<<"In Total We Have "<<connectedMap.size()<<" Parts"<<std::endl;
			//Mask=NULL;

		}


	

}