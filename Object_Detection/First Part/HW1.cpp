//Read & WRITE  Video from file using VideoCapture and VideoWriter and GRAB and RETRIEVE
//getPerspective and warpPerspective with image and video

#include<opencv/cv.h>
#include<opencv/highgui.h>
#include <iostream>
#include <opencv2/opencv.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>
#include <locale>
#include <cctype>

using namespace	 std;
using namespace cv;

vector<Point2f> inputPoints;
vector<Point2f> outputPoints;
int WRONG_POINTS=0;


void parseInputTextFile(char* argv[]){

	  //text file
	  std::string file;
	  file=argv[2];

	   //create a file-reading object
	  ifstream fin;
	  fin.open(file.c_str(),std::ifstream::in);
	  int flag=0;
	  //Check if file is open
	  if (!fin.good())
	   cout<<"Could not open text file\n";
	  else{
		  std::string str;
    		//Read file and save the points

		  while (!fin.eof())
		   {
			  getline(fin,str);

			 if(flag!=4){

		 	 std::stringstream ss;
		 	 ss<<str;
		 	 float i1=0,i2=0,o1=0,o2=0;
		 	 ss>>i1>>i2>>o1>>o2;

		 	 inputPoints.push_back(Point(i1,i2));
		 	 outputPoints.push_back(Point(o1,o2));

		 	 flag++;
			 }
		   }
		  if(flag==3)
		   {
		       	 cout<<"need 4 points to generate perspective transformation";
		       	 WRONG_POINTS=1;
		   }

	  }
}

int perspectTransform(char* argv[]){

	//initialize capture to read input video file
	VideoCapture cap;
	cap.open(argv[1]);

	if(!cap.isOpened())  // check if open
        return -1;

	//create window to show video
	namedWindow("Original",1);
	double w = cap.get(CV_CAP_PROP_FRAME_WIDTH); //get the width of frames of the video
    double h = cap.get(CV_CAP_PROP_FRAME_HEIGHT); //get the height of frames of the video

	//get the frame size
   Size size(static_cast<int>(w), static_cast<int>(h));

   //Give the output filename to the video

   std::stringstream outfile;
   outfile << argv[3] << ".mpg";
   std::string fileName = outfile.str();


   //initialize the VideoWriter object to save the output of the tranformed video
   VideoWriter out (fileName, CV_FOURCC('M','P','E','G'), cap.get(CV_CAP_PROP_FPS), size, true);

    if ( !out.isOpened() ) //if not initialize the outWriter successfully, exit the program
   {
		cout << "ERROR: Failed to write the video" << endl;
		return -1;
   }


   //Create matrix to store frame
   Mat frame;

   while(1){

		bool ifgrab=cap.grab();
		if(!ifgrab)
			break;
		bool ifwrite=cap.retrieve(frame);
		if(!ifwrite)
			break;

		imshow("Original",frame);

		//getPerspective and warpPerspective with video

		//cout<<"Input"<<endl;
		//cout<<inputPoints<<endl;
		//cout<<"output"<<endl;

	    //Computes Perspective transformation
		Mat trans=getPerspectiveTransform(inputPoints,outputPoints);

		//Transform the image
		Mat transform(frame.size(),frame.type());
		warpPerspective(frame,transform,trans,Size(w,h));

		//Write the frame to videowriter
		out<<transform;

		//print frame to screen
		imshow("Transformed Video",transform);

		//delay 30ms
		waitKey(30);
	}

	//Release the objects
	out.release();
	cap.release();

	return 0;
}

int main(int argc, char* argv[]){

	//Check if exactly 3 Arguments are provided
	if( argc != 4)
    {
	 cout<<"Requires three (name and path) input: TextFile Video OutputFilename (without extension)"<<endl;
     return -1;
    }

	//Process the text file and save the points
	parseInputTextFile(argv);
	if(WRONG_POINTS==1)
				return -1;
	//Process the input video and generate output
	int flag=perspectTransform(argv);
	if(flag==-1)
		return -1;

	return 0;

}
