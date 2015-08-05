#include<opencv/cv.h>
#include<opencv/highgui.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include<fstream>

using namespace std;
using namespace cv;

void findMatch(const Mat& input, vector<Point>& loc, size_t size) {
	float maxValue = -1.0f * numeric_limits<float>::max();
	float* data = reinterpret_cast<float*>(input.data);

	for (int i = 0; i < input.rows; i++) {
		for (int j = 0; j < input.cols; j++) {

			if (data[i * input.cols + j] > maxValue) {
				maxValue = data[i * input.cols + j];

				if (loc.size() < (size + 1)) {
					loc.push_back(Point(j, i));
				}
			}
		}
	}
}

int main(int argc, char* argv[]) {

	//Check if exactly 2 Arguments are provided
	if (argc != 3) {
		cout
				<< "Requires three (name and path) input: TextFile Video OutputFilename (without extension)"
				<< endl;
		return -1;
	}
	cv::Mat template_img;
	cv::Mat result_mat2;
	cv::Mat debug_img2;

	template_img = cv::imread(argv[2], CV_LOAD_IMAGE_GRAYSCALE);
	if (template_img.data == NULL) {
		printf("imread() failed...\n");
		return -1;
	}

	//initialize capture to read input video file
	VideoCapture cap;
	cap.open(argv[1]);

	if (!cap.isOpened())  // check if open
		return -1;

	//create window to show video

    namedWindow("RESULT", CV_WINDOW_AUTOSIZE);
	double w = cap.get(CV_CAP_PROP_FRAME_WIDTH); //get the width of frames of the video
	double h = cap.get(CV_CAP_PROP_FRAME_HEIGHT); //get the height of frames of the video
	//cout << "Frame Size = " << w << "x" << h << endl;

	//get the frame size
	Size size(static_cast<int>(w), static_cast<int>(h));


	//Output File
	ofstream ofile ("objectDetected.txt");
	if(!ofile.is_open())
	{
	    cout<<"Could not open the file"<<endl;
	    return -1;
	}

	Mat frame; //Declare frame to grab video
	int frame_no = 0; //Declare variable to keep track of frame number

	while (1) {

		bool ifgrab = cap.grab();
		if (!ifgrab)
			break;
		bool ifwrite = cap.retrieve(frame);
		if (!ifwrite)
			break;
		cv::cvtColor(frame, debug_img2, CV_BGR2GRAY);



		int match_method = CV_TM_SQDIFF_NORMED; //CV_TM_CCORR_NORMED;
		cv::matchTemplate(debug_img2, template_img, result_mat2, match_method);

		//threshold(result_mat2,result_mat2,0.9,1.,CV_THRESH_TOZERO);
		//normalize(result_mat2, result_mat2, 0, 1, NORM_MINMAX, -1, Mat());

		double minVal;
		double maxVal;
		Point minLoc, maxLoc, matchLoc;

		//////////////////////////////////Handle Multiple Occurrences///////////////////////////////////////////////////
		Mat old_mat = result_mat2.clone();
		//Image Pyramid
		// Mat imgPyr;
		//pyrUp( template_img, imgPyr, Size( template_img.cols*2, template_img.rows*2 ) );
		//imshow("PYR",imgPyr);
		/*
		/// For SQDIFF and SQDIFF_NORMED, the best matches are lower values. For all the other methods, the higher the better
		if (match_method == CV_TM_SQDIFF
				|| match_method == CV_TM_SQDIFF_NORMED) {
			result_mat2 = 1.0 - result_mat2;
		}
		// get the top 10 maximums...
		vector<Point> res;
		findMatch(result_mat2, res, 10);

		int i = 1;
		while (i < res.size()) {
			Point matchLoc = res.at(i);

			rectangle(frame, matchLoc,
						Point(matchLoc.x + template_img.cols,
								matchLoc.y + template_img.rows),
						CV_RGB(0, 0, 255), 3);
			ofile<<"Swiss Ball"<<frame_no<<" "<<matchLoc.x<<" "<<matchLoc.y<<" "<<(matchLoc.x + template_img.cols)<<" "<<(matchLoc.y + template_img.rows);
			ofile<<"\n";

			i++;

		}
		*/
		////////////////////////////////////////////Threshold the value and display//////////////////////////////////////////

		minMaxLoc(old_mat, &minVal, &maxVal, &minLoc, &maxLoc, Mat());
		//cout<<minVal<<" "<<maxVal<<endl;
		if (match_method == CV_TM_SQDIFF || match_method == CV_TM_SQDIFF_NORMED)
			matchLoc = minLoc;
		else
			matchLoc = maxLoc;

		if(minVal<0.15 && minVal>0)	{

			rectangle(frame, matchLoc,
					Point(matchLoc.x + template_img.cols,
							matchLoc.y + template_img.rows), CV_RGB(0, 0, 255),
					3);

		}

        imshow("RESULT", frame);
		//Object frame UpperLeftX UpperLeftY LowerRightX LowerRightY
		ofile<<"Swiss Ball "<<frame_no<<" "<<matchLoc.x<<" "<<matchLoc.y<<" "<<(matchLoc.x + template_img.cols)<<" "<<(matchLoc.y + template_img.rows);
		ofile<<"\n";
		frame_no++;

		//delay 30ms
		waitKey(30);
	}
	//Release the objects
	cap.release();

	return 0;
}
