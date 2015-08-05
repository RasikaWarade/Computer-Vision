#include<opencv\cv.h>
#include<opencv\highgui.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/nonfree/features2d.hpp"
#include "opencv2/nonfree/nonfree.hpp"
#include <math.h>
#include <fstream>

using namespace std;
using namespace cv;

Mat full_dictionary;
int no_videos=7;//*** change the number of videos as per train data set
int N=10;//number of frames ***

int imagesData=0;
float idfVector[100][1000];

void idf_vector(Mat full_dictionary){
	ofstream myfile;
	myfile.open ("example.txt");
	myfile << "Calculating IDF_VECTORS.\n";
	
	std:: string videoName="";
	
	int n_frames[100];
	//create dictionary
	int dict_size=100;//***
	

	//create a nearest neighbor matcher
    Ptr<DescriptorMatcher> matcher(new FlannBasedMatcher);
    //create Sift feature point extracter
    Ptr<FeatureDetector> detector(new SiftFeatureDetector());
    //create Sift descriptor extractor
    Ptr<DescriptorExtractor> extractor(new SiftDescriptorExtractor);    
    //create BoF (or BoW) descriptor extractor
    BOWImgDescriptorExtractor bowDE(extractor,matcher);
    //Set the dictionary with the vocabulary we created in the first step
    bowDE.setVocabulary(full_dictionary);
	
	for(int i=1;i<no_videos;i++){ 
		
		stringstream temp;
		temp<<i;
		std::string no=temp.str();
		videoName="C:/Rasika/video_"+no+".avi"; //*** path can be changed 
		
		//initialize capture
		VideoCapture cap;
		cap.open(videoName);
		
		double count = cap.get(CV_CAP_PROP_FRAME_COUNT); //get the frame count
		
		int jump=count/N; //extract 10 frames from the video ***
		int j=0;
		int cnt=0;
		myfile<<"Reading Video";
		Mat features;
		Mat desc;
		while(cnt<count){
			
			//Create matrix to store video frame
			Mat image;
			cap.set(CV_CAP_PROP_POS_FRAMES,cnt); //Set index to jump for particular count
			bool success = cap.read(image); 
			if (!success){
				cout << "Cannot read  frame " << endl;
				break;
			}
			
			///////////Convert to gray scale/////////////
			Mat gray_image;
			cvtColor( image, gray_image, CV_BGR2GRAY );
			imagesData++;//Number of images in the database

			 //To store the keypoints that will be extracted by SIFT
			vector<KeyPoint> keypoints;        
			//Detect SIFT keypoints (or feature points)
			detector->detect(gray_image,keypoints);
			//To store the BoW (or BoF) representation of the image
			Mat bowDescriptor;        
			//extract BoW (or BoF) descriptor from given image
			bowDE.compute(gray_image,keypoints,bowDescriptor);
			
			desc.push_back(bowDescriptor);
					
			////////////////
			//delay 33ms //***
			//waitKey(33);
		
			cnt+=jump;
			j++;

			///next frame for the same video
		}

		
		
			/*myfile<<desc.rows<<endl;
			myfile<<desc.cols<<endl;

			int tf=0;
			for(int i=0;i<desc.rows;i++){
				for(int j=0;j<desc.cols;j++){
					if(desc.at<float>(i,j)>0){

						//cout<<bowDescriptor.at<float>(i,j)<<endl;
						tf++;
					}
				}
			}
	
			myfile<<"Term Frequency:"<<tf<<"\n";
			float idf=0;
			float logcal=count/tf;
			idf=log(logcal);
			myfile<<"IDF:"<<idf<<"\n";
			idfVector[i-1][j]=idf;
			
			myfile<<idfVector[i-1][j];*/
			
		//store number of frames per video
		n_frames[i-1]=j;
	
		
	}
	myfile<<"IDF done";
	myfile.close();

	

}
int trainData(){
		
	std:: string videoName="";
	
	int n_frames[1000];
	//create dictionary
	int dict_size=100;//***
	
	Mat features;
	for(int i=1;i<no_videos;i++){ 
		
		
		stringstream temp;
		temp<<i;
		std::string no=temp.str();
		videoName="C:/Rasika/trainvideos/video_"+no+".avi"; //*** path can be changed 
		
		//initialize capture
		VideoCapture cap;
		cap.open(videoName);
		if(!cap.isOpened())  // check if we succeeded
			return -1;
		
		double count = cap.get(CV_CAP_PROP_FRAME_COUNT); //get the frame count
		
		//create window to show image
		//namedWindow("Video",1);
		//cout<<count<<endl;
		int jump=count/N; 
		int j=1;
		
		int u=0;
		if(count<10){
			jump=1;
		}
		int cnt=jump;
		while(u<10){
			
			//Create matrix to store video frame
			Mat image;
			cap.set(CV_CAP_PROP_POS_FRAMES,cnt); //Set index to jump for particular count
			bool success = cap.read(image); 
			if (!success){
				cout << "Cannot read  frame " << endl;
				break;
			}
			
			///////////Convert to gray scale/////////////
			Mat gray_image;
			cvtColor( image, gray_image, CV_BGR2GRAY );

			////////EXTRACT INTEREST POINTS USING SIFT////
			// vector of keypoints
			std::vector<cv::KeyPoint> keypoints;
			// Construct the SIFT feature detector object
			SiftFeatureDetector sif(0.03,10.); // threshold  //***
			//Detect interest points
			sif.detect(gray_image,keypoints);

			////////IMSHOW THE FRAMES EXTRACTED///////////
			
			//copy video stream to image
			//cap>>image;
			//print image to screen
			//imshow("Video",image);
			

			///////////Save the frames//////////////
			
			stringstream temp2;
			temp2<<j;
			std::string no2=temp2.str();
			std::string frame_name="frame"+no2+".jpg";
			imwrite(frame_name,image);
			
			
			//////////////Draw the keypoints////////////
			
			/*
			Mat featureImage;
			// Draw the keypoints with scale and orientation information
			drawKeypoints(image, // original image
			keypoints, // vector of keypoints
			featureImage, // the resulting image
			Scalar(255,0,255), // color of the points
			DrawMatchesFlags::DRAW_RICH_KEYPOINTS); //flag
			//std::string name="image"+i;
			imshow(frame_name, featureImage );
			*/
			
			////////////////////detect decriptors//////////////////
			
			SiftDescriptorExtractor siftExtractor;
			Mat siftDesc;
			siftExtractor.compute(gray_image,keypoints,siftDesc);
			features.push_back(siftDesc);//add the descriptors from each frame..to create one for a video
		
			////////////////
			//delay 33ms //***
			//waitKey(33);
			
			cnt+=jump;
			j++;
			u++;
			///next frame for the same video
		}
			
		//store number of frames per video
		n_frames[i-1]=j-1;
	
			
		
	}
	
		TermCriteria term(CV_TERMCRIT_ITER,100,0.001);//***
		
		//retries number ***
		int retries=1; 
		
		int flags=KMEANS_PP_CENTERS;
		BOWKMeansTrainer bowTrainer(dict_size,term,retries,flags);
		//cluster the feature vectors
		Mat dictionary=bowTrainer.cluster(features); 

		//for further process
		full_dictionary.push_back(dictionary);
		///////////////////////////////////////////////////
		FileStorage fs("full_dictionary.yml", FileStorage::WRITE);
		fs << "vocabulary" << full_dictionary;
		fs.release();
		//Created Vocabulary

		//Calculate histograms for the train videos
		//idf_vector(full_dictionary);
		
		return 0;
}

int main(int argc, char** argv){

	//Check arguments
	//***
	
	int set=trainData();
	 
	//If set=0, proceed
	if(set==0){
	

	//Take the two video inputs and measure the similarity
	float firstTF[1000];//***
	float secondTF[1000];
	int n_frames[1000];

	//////////////////////////////////////////////////////////////////////////////////////////////////	
	Mat dicty; 
	FileStorage fs("full_dictionary.yml", FileStorage::READ);
    fs["vocabulary"] >> dicty;
    fs.release();

	//set dictionary
	int dict_size=100;//***
	//create a nearest neighbor matcher
    Ptr<DescriptorMatcher> matcher(new FlannBasedMatcher);
    //create Sift feature point extracter
    Ptr<FeatureDetector> detector(new SiftFeatureDetector());
    //create Sift descriptor extractor
    Ptr<DescriptorExtractor> extractor(new SiftDescriptorExtractor);    
    //create BoF (or BoW) descriptor extractor
    BOWImgDescriptorExtractor bowDE(extractor,matcher);
    //Set the dictionary with the vocabulary we created in the first step
    bowDE.setVocabulary(dicty);
	
//////////////////////////////First Video//////////////////////////////////////////////////////////
	
	ofstream myfile;
	myfile.open ("first_video.txt");
	myfile << "Calculating TF_VECTORS.\n";

		//initialize capture
		VideoCapture cap;
		cap.open(argv[1]); //***
		
		double count = cap.get(CV_CAP_PROP_FRAME_COUNT); //get the frame count
		
		int jump=count/N; //extract 10 frames from the video ***
		int j=0;
		if(count<10){
			jump=1;
		}
		int cnt=jump;
		myfile<<"Reading Video";
		Mat features;
		Mat desc;

		int u=0;
		while(u<10){
			
			//Create matrix to store video frame
			Mat image;
			cap.set(CV_CAP_PROP_POS_FRAMES,cnt); //Set index to jump for particular count
			bool success = cap.read(image); 
			if (!success){
				cout << "Cannot read  frame " << endl;
				break;
			}
			
			///////////Convert to gray scale/////////////
			Mat gray_image;
			cvtColor( image, gray_image, CV_BGR2GRAY );

			 //To store the keypoints that will be extracted by SIFT
			vector<KeyPoint> keypoints;        
			//Detect SIFT keypoints (or feature points)
			detector->detect(gray_image,keypoints);
			//To store the BoW (or BoF) representation of the image
			Mat bowDescriptor;        
			//extract BoW (or BoF) descriptor from given image
			bowDE.compute(gray_image,keypoints,bowDescriptor);
			
			desc.push_back(bowDescriptor);
			
			cnt+=jump;
			j++;
			u++;
			///next frame for the same video
		}
		//FileStorage fs("descriptor.yml", FileStorage::WRITE);
		//fs << "descriptor" << desc;
		//fs.release();
		
		
		
		for(int k=0;k<desc.cols;k++){
			int tf=0;
			for(int l=0;l<desc.rows;l++){
					if(desc.at<float>(l,k)>0){

						//cout<<bowDescriptor.at<float>(i,j)<<endl;
						tf++;
					}


				}
			myfile<<"Term Frequency:"<<tf<<"\n";
			firstTF[k]=tf;
			
		}
	
		
			
		myfile<<"TF done";
		myfile.close();

//////////////////////////////Second Video//////////////////////////////////////////////////////////
		
		ofstream myfile3;
		myfile3.open ("second_video.txt");
		myfile3 << "Calculating IDF_VECTORS.\n";

		//initialize capture
		cap.open(argv[2]); //***
		
		count = cap.get(CV_CAP_PROP_FRAME_COUNT); //get the frame count
		
		jump=count/N; //extract 10 frames from the video ***
		j=0;
		if(count<10){
			jump=1;
		}
		cnt=jump;
		myfile3<<"Reading Video";
		Mat desc2;
		u=0;
		while(u<10){
			
			//Create matrix to store video frame
			Mat image;
			cap.set(CV_CAP_PROP_POS_FRAMES,cnt); //Set index to jump for particular count
			bool success = cap.read(image); 
			if (!success){
				cout << "Cannot read  frame " << endl;
				break;
			}
			
			///////////Convert to gray scale/////////////
			Mat gray_image;
			cvtColor( image, gray_image, CV_BGR2GRAY );

			 //To store the keypoints that will be extracted by SIFT
			vector<KeyPoint> keypoints;        
			//Detect SIFT keypoints (or feature points)
			detector->detect(gray_image,keypoints);
			//To store the BoW (or BoF) representation of the image
			Mat bowDescriptor;        
			//extract BoW (or BoF) descriptor from given image
			bowDE.compute(gray_image,keypoints,bowDescriptor);
			
			desc2.push_back(bowDescriptor);
			cnt+=jump;
			j++;
			u++;
			///next frame for the same video
		}
			

		for(int k=0;k<desc2.cols;k++){
			int tf=0;
			for(int l=0;l<desc2.rows;l++){
					if(desc2.at<float>(l,k)>0){

						//cout<<bowDescriptor.at<float>(i,j)<<endl;
						tf++;
					}


				}
			myfile3<<"Term Frequency:"<<tf<<"\n";
			secondTF[k]=tf;
			
		}
			
			myfile3<<"TF done";
			myfile3.close();

//////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//Display the similarity score 
		
		//Dot product of TF vectors 

		
		float similarity=0;
		ofstream my3;
		my3.open("Similarity.txt");
		
		for(int i=0;i<dict_size;i++){
				similarity+=firstTF[i]*secondTF[i];
				
		}
			my3<<"\n";
		my3<<similarity<<" ";
		my3.close();

		cout<<"Similarity Score:"<<similarity<<endl;
		
	}


	return 0;
}