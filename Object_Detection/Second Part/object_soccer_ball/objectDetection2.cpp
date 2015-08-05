#include<opencv/cv.h>
#include<opencv/highgui.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include<fstream>
#include "opencv2/imgproc/imgproc.hpp"

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
void PCA_method(char* argv[]) {

	Mat template1 = imread("temp1.JPG", CV_LOAD_IMAGE_GRAYSCALE); //read the image data in the file "MyPic.JPG" and store it in 'img'
	Mat template2 = imread("temp2.JPG", CV_LOAD_IMAGE_GRAYSCALE);
	Mat template3 = imread("temp3.JPG", CV_LOAD_IMAGE_GRAYSCALE);
	Mat main_template = imread("template_42.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	if (template1.empty()) //check whether the image is loaded or not
	{
		cout << "Error : Image cannot be loaded..!!" << endl;
		//system("pause"); //wait for a key press

	}

	//Mat first1=main_template(Rect(268,36,38,38)).clone();
	//imshow("actual",first1);
	normalize(template1, template1, 0, 1, NORM_MINMAX, CV_64F); //normalize. the intensity value is between 0 and 1 for all pixels
	normalize(template2, template2, 0, 1, NORM_MINMAX, CV_64F);
	normalize(template3, template3, 0, 1, NORM_MINMAX, CV_64F);
	normalize(main_template, main_template, 0, 1, NORM_MINMAX, CV_64F);

	//template2=template2;
	//template3=template3;
	//cout<<template1;

	///////////////////////////////////////////// RESHAPE THE TEMPLATE IMAGES TO COLUMN VECTORS ///////////////////////////

	Mat img1_template1 = template1.reshape(0, 1);
	Mat img2_template2 = template2.reshape(0, 1);
	Mat img3_template3 = template3.reshape(0, 1);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Mat img3_template3=template1.reshape(0,template3.rows*template3.cols);

	///////////////////////////////////////////// CALCULATE MEAN OF EACH TEMPLATE AND SUBTRACT /////////////////////////////////////////////
	Scalar tempVal1 = mean(template1);
	Scalar tempVal2 = mean(template2);
	Scalar tempVal3 = mean(template3);

	double myMatMean1 = tempVal1.val[0];
	double myMatMean2 = tempVal2.val[0];
	double myMatMean3 = tempVal3.val[0];

	cout << template1.rows << '\t' << template1.cols << endl;
	cout << template2.rows << '\t' << template2.cols << endl;
	cout << template3.rows << '\t' << template3.cols << endl;

	cout << "Mean 1: " << myMatMean1 << endl;
	cout << "Mean 2: " << myMatMean2 << endl;
	cout << "Mean 3: " << myMatMean3 << endl;

	Mat resultant1, resultant2, resultant3;

	subtract(myMatMean1, img1_template1, resultant1);
	subtract(myMatMean2, img2_template2, resultant2);
	subtract(myMatMean3, img3_template3, resultant3);

	cout << "Dimensions of the three templates as column matrices" << endl;
	cout << resultant1.rows << '\t' << resultant1.cols << endl;
	cout << resultant2.rows << '\t' << resultant2.cols << endl;
	cout << resultant3.rows << '\t' << resultant3.cols << endl;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////// CONCATENATE COLUMN MATRICES AS ONE ////////////////////////////////////////////////

	/*
	 Mat final_template=resultant1.clone();
	 final_template.push_back(resultant2);
	 final_template.push_back(resultant3);

	 final_template=final_template.reshape(0,1444);
	 cout<<final_template.rows<<'\t'<<final_template.cols<<endl;
	 */

	///*
	Mat final_template = img1_template1.clone();
	final_template.push_back(img2_template2);
	final_template.push_back(img3_template3);

	final_template = final_template.reshape(0, 1444);
	cout << final_template.rows << '\t' << final_template.cols << endl;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////// CALCULATE MEAN OF FINAL TEMPLATE AND SUBTRACT /////////////////////////////////////////////

	Scalar tempVal11 = mean(final_template);
	Mat final_resultant;
	double myMatMean_final = tempVal11.val[0];
	subtract(myMatMean_final, final_template, final_resultant);

	//cout<<final_resultant<<endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////// CALCULATE PCA FOR EIGEN VALUES AND EIGEN VECTORS /////////////////////////////////////////////

	PCA pca(final_resultant, // pass the data
			Mat(), // there is no pre-computed mean vector,
				   // so let the PCA engine to compute it
			CV_PCA_DATA_AS_COL, // indicate that the vectors
								// are stored as matrix rows
								// (use CV_PCA_DATA_AS_COL if the vectors are
								// the matrix columns)
			2 // specify how many principal components to retain
			);

	Mat resultant4;
	Mat eigenvalues = pca.eigenvalues.clone();
	Mat eigenvectors = pca.eigenvectors.clone();
	Mat mean_val= pca.mean.clone();
	//subtract(mean_val, main_template, resultant4);

	Mat output1;
	cout << endl << endl << "The dimension of eigen vectors " << endl;
	cout << eigenvectors.rows << '\t' << eigenvectors.cols << endl;
	cout << eigenvalues << endl; //'\t'<<mean_val<<endl;
	//pca.project(resultant4, output1);

	///////////////////////////////////////// TRANSPOSE OF THE EIGEN VECTOR MATRIX /////////////////////////////////////////

	Mat new_img;
	transpose(eigenvectors, new_img);
	//cout<<new_img.rows<<'\t'<<new_img.cols<<endl;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////// RESHAPING EIGEN VECTOR IN COLUMN 1 ////////////////////////////////////////////////

	Mat eigen_first = eigenvectors.row(0);
	Mat eigen_second = eigenvectors.row(1);
	Mat eigen_third = eigenvectors.row(1);
	eigen_first = eigen_first.reshape(0, 38);
	eigen_second = eigen_second.reshape(0, 38);
	eigen_third = eigen_third.reshape(0, 38);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////// FOR MAIN IMAGE /////////////////////////////////////////////////////////////

	/////////////////////////////// main_template: Normalized image //////////////////////////////////////////////////
	Mat main_image;
	//main_image=main_template-myMatMean_final;
	Scalar tempVal12 = mean(main_template);
	double final_image_resultant = tempVal12.val[0];

	subtract(myMatMean_final, main_template, main_image);

	//main_image=main_image.reshape(0,2);
	//transpose(main_image,main_image);

	mean_val = mean_val.reshape(0, 38);
	cout << mean_val.rows << '\t' << mean_val.cols << endl;
	imshow("mean", mean_val);

	Mat result;
	//matchTemplate(main_image, eigen_first, result, CV_TM_CCOEFF_NORMED );

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Rect rect;
	Mat first = main_image(Rect(268, 36, 38, 38)).clone();
	cout << eigen_first << endl;
	double dist;
	double dist_2;
	dist = first.dot(eigen_first);
	cout << "dist" << dist << endl;
	dist_2 = first.dot(eigen_second);
	cout << "dist" << dist_2 << endl;

	Mat sec = main_image(Rect(290, 80, 38, 38)).clone();
	//imshow("!",sec);
	double dist1;
	dist1 = sec.dot(eigen_first);
	cout << "dist" << dist1 << endl;
	double dist1_2;
	dist1_2 = sec.dot(eigen_second);
	cout << "dist" << dist1_2 << endl;
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << "First eigen value:" << eigenvalues.row(0);

	double val1 = eigenvalues.at<double>(0, 0);
	cout << val1 << endl;
	double product1 = 1;
	double product2 = 1;

	double product1_2 = 1;
	double product2_2 = 1;

	//First eigenvalue and two eigenvectors//
	double rho = sqrt(val1);
	double up = -(dist * dist) / (2 * val1);
	double res = (1 / (rho * sqrt(2 * 3.14))) * exp(up);
	cout << res << endl;
	product1 = product1 * res;

	up = (dist_2 * dist_2) / (2 * val1);
	res = (1 / (rho * sqrt(2 * 3.14))) * exp(up);

	product1 = product1 * res;

	up = -(dist1 * dist1) / (2 * val1);
	res = (1 / (rho * sqrt(2 * 3.14))) * exp(up);
	cout << res << endl << endl;
	product2 = product2 * res;

	up = -(dist1_2 * dist1_2) / (2 * val1);
	res = (1 / (rho * sqrt(2 * 3.14))) * exp(up);
	cout << res << endl << endl;
	product2 = product2 * res;

	eigen_first = eigen_first * val1;
	imshow("image", eigen_first);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/*
	 dist=first.dot(eigen_second);
	 //cout<<"dist"<<dist<<endl;


	 dist1=sec.dot(eigen_second);
	 //float dd=ceil(dist1*100)/100;
	 imshow("abc",first);
	 imshow("second",sec);
	 //cout<<"dist"<<dist1<<endl;

	 //cout<<"First eigen value:"<<eigenvalues.row(0);

	 val1 = eigenvalues.at<double>(1,0);
	 //cout<<val1<<endl;


	 rho=sqrt(val1);
	 up=-(dist*dist)/(2*val1);
	 res=(1/(rho*sqrt(2*3.14)))*exp(up);
	 cout<<res<<endl;
	 product1=product1*res;
	 up=-(dist1*dist1)/(2*val1);
	 res=(1/(rho*sqrt(2*3.14)))*exp(up);
	 cout<<res;
	 product2=product2*res;

	 */
	cout << "Prob of ball: " << product1 << endl;
	cout << "Prob of no ball: " << product2 << endl;

	//////////////////////////////////////////////////////////////////////
	//imshow("source_window",main);
	Mat result_mat2;
	int match_method = CV_TM_CCORR_NORMED;/*
	 cv::matchTemplate(main_template, eigen_first, result_mat2, match_method);

	 cv::normalize(result_mat2, result_mat2, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());

	 double minVal; double maxVal;
	 cv::Point minLoc, maxLoc, matchLoc;
	 cv::minMaxLoc(result_mat2, &minVal, &maxVal, &minLoc, &maxLoc, cv::Mat() );
	 if( match_method == CV_TM_SQDIFF || match_method == CV_TM_SQDIFF_NORMED ) matchLoc = minLoc;
	 else matchLoc = maxLoc;
	 cvtColor(main_template, main_template, CV_GRAY2BGR);
	 cv::rectangle(
	 main_template,
	 matchLoc,
	 cv::Point(matchLoc.x + eigen_first.cols , matchLoc.y + eigen_first.rows),
	 CV_RGB(255,0,0),
	 3);
	 imshow( "result_window", main_template );

	 */

	/*
	 namedWindow("template", CV_WINDOW_AUTOSIZE);
	 namedWindow("template1", CV_WINDOW_AUTOSIZE);
	 namedWindow("template2", CV_WINDOW_AUTOSIZE);
	 namedWindow("template3", CV_WINDOW_AUTOSIZE);//create a window with the name "MyWindow"
	 imshow("template1", template1); //display the image which is stored in the 'img' in the "MyWindow" window
	 imshow("template2", template2);
	 //imshow("template3", result);
	 imshow("template", main_image);
	 */

	imshow("template", mean_val);
	waitKey(0); //wait infinite time for a keypress

	destroyWindow("MyWindow"); //destroy the window with the name, "MyWindow"


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

	///////////////////////////PCA/////////////////////////
	//Provide more than one templates to run PCA
	//PCA_method(argv);
	///////////////////////////////////////////////////////

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
	ofstream ofile("objectDetected.txt");
	if (!ofile.is_open()) {
		cout << "Could not open the file" << endl;
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
		 ofile<<"Soccer Ball"<<frame_no<<" "<<matchLoc.x<<" "<<matchLoc.y<<" "<<(matchLoc.x + template_img.cols)<<" "<<(matchLoc.y + template_img.rows);
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

		if (minVal < 0.07) {//0.07 if not zero mean and normalized//normalized-minVal<0.7 && minVal>0.1

			rectangle(frame, matchLoc,
					Point(matchLoc.x + template_img.cols,
							matchLoc.y + template_img.rows), CV_RGB(0, 0, 255),
					3);

		}

		imshow("RESULT", frame);
		//Object frame UpperLeftX UpperLeftY LowerRightX LowerRightY
		ofile << "Soccer Ball " << frame_no << " " << matchLoc.x << " "
				<< matchLoc.y << " " << (matchLoc.x + template_img.cols) << " "
				<< (matchLoc.y + template_img.rows);
		ofile << "\n";
		frame_no++;

		//delay 30ms
		waitKey(30);
	}
	//Release the objects
	cap.release();

	return 0;
}
