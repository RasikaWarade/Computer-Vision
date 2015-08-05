//////////////////////////////////////////////////////Headers///////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>
#include <locale>
#include <cctype>
#include <vector>

/////////////////////////////////////////////////////Declare Global variable////////////////////////////////////////////////////////////////////////////////////

using namespace std;
int iFrames, oFrames, frameLinks; //Numer of frames in Input and output frames and total number of links for robot
std::vector<std::string> words;
double position[1000][1000], velocity[1000][1000]; //Arrays to store position and velocity of control points of robot

////////////////////////////////////////////////////Define Functions///////////////////////////////////////////////////////////////////////////////////

/*
 * Func: Parse the input text file and store the values in the variables
 */

void parseInputTextFile(char* argv[]) {

	bool readFrames = true;
	bool readLinks = true;

	//text file
	std::string file = "robot.key"; //***
	//file = argv[1]; ***

	//create a file-reading object
	ifstream fin;
	fin.open(file.c_str(), std::ifstream::in);

	//Check if file is open
	if (!fin.good())
		cout << "Could not open robot.key file\n";
	else {
		std::string str;

		//Read file and save the points
		int p_cnt = 0;
		while (!fin.eof()) {
			getline(fin, str);

			std::stringstream ss;
			ss << str;
			if (str.length() != 0) {

				if (readFrames) {
					ss >> iFrames >> oFrames;
					if (iFrames != 0 && oFrames != 0) {
						readFrames = false;
					}
				} else {

					if (readLinks) {

						//Check no. of links
						//take first frame details
						double f;
						double pt[1000];
						int i = 0;
						std::string word;

						while (ss >> f) {
							pt[i] = f;
							word = "1";
							words.push_back(word);
							i++;
						}
						frameLinks = words.size();

						if (iFrames != 0 && frameLinks != 0) {
							for (int j = 0; j < frameLinks; j++) {
								position[p_cnt][j] = pt[j];
							}

							getline(fin, str);
							std::stringstream ss2;
							ss2 << str;

							double v;
							int k = 0;
							while (ss2 >> v) {
								velocity[p_cnt][k] = v;
								k++;
							}
							readLinks = false;
							p_cnt++;
						}

					} else {

						//take rest frames
						double v;
						int t = 0;
						while (ss >> v) {
							position[p_cnt][t] = v;
							t++;
						}

						getline(fin, str);
						std::stringstream ss2;
						ss2 << str;
						int y = 0;
						while (ss2 >> v) {
							velocity[p_cnt][y] = v;
							y++;
						}

						p_cnt++;

					}

				}
			}
		}

	}
}

/*
 * Func: Find the P(u) using Hermite Equation
 */

double findQ(double s[1][4], double H[4][4], double C[4][1]) {
	double first[4];
	double t = 0;
	///////Matrix Multiplication////////////
	//Find P(u)
	// 1x4 multiply 4x4
	for (int r = 0; r < 4; r++) {
		t = 0;
		for (int c = 0; c < 4; c++) {

			t = t + s[0][c] * H[c][r];

		}

		first[r] = t;
	}

	//P(u)=q
	double q = 0;

	//1x4 multiply 4x1
	for (int i = 0; i < 4; i++) {
		q = q + first[i] * C[i][0];
	}
	////////////////////////////////////////
	return q;
}

/*
 * Func: Perform Position Interpolation on the given control points
 */
void positionInterpolation() {

	double interpolate[100][100];

	//calculate delta
	double delta = (double) ((iFrames - 1)) / (oFrames - 1);
	int noCurves = iFrames - 1;
	int framesMinus = oFrames - iFrames;
	int forEachCurve = (framesMinus - 2) / noCurves;

	//Create hermite matrix
	double H[4][4] = { { 2, -2, 1, 1 }, // row 0
			{ -3, 3, -2, -1 }, // row 1
			{ 0, 0, 1, 0 }, // row 2
			{ 1, 0, 0, 0 } //row 3
	};

	double s[1][4] = { 1, 1, 1, 1 };
	double C[4][1] = { { 1 }, { 1 }, { 1 }, { 1 } };

	//for each control point, for each frame interpolate results
	//initializing C and s
	//calulate interpolate frame at rest and add frame 0

	//Assume it starts and ends at rest --Total 2 Frames

	for (int i = 0; i < frameLinks; i++) {
		C[0][0] = position[0][i];
		C[1][0] = position[0][i];
		C[2][0] = 0;
		C[3][0] = velocity[0][i];
		double q = 0;
		double s_delta = delta;
		s_delta = 1;
		s[0][1] = s_delta * s_delta * s_delta;
		s[0][2] = s_delta * s_delta;
		s[0][3] = s_delta;
		s[0][4] = 1;
		q = findQ(s, H, C);
		interpolate[0][i] = q;

	}


	//include first frame to output--Total 1 Frame
	for (int p = 0; p < frameLinks; p++) {
		interpolate[1][p] = position[0][p];
	}

	int of = 2;	//start for output second frame
	int delta_count;

	for (int i = 0; i < frameLinks; i++) {	//--Total 64 frames
		//take two frames calculate C
		//calculate s then S[][]

		delta_count = 2;
		of = 2;
		for (int j = 0; j < iFrames - 1; j++) {	//5 frames then 4 curves so repeat 4 times

			//Calculate C
			C[0][0] = position[j][i];
			C[1][0] = position[j + 1][i];
			C[2][0] = velocity[j][i];
			C[3][0] = velocity[j + 1][i];

			//for each frame calculate intermediate 15 frames and add last one
			for (int l = 0; l < forEachCurve; l++) {

				double s_delta = 0;
				s_delta = delta * delta_count;
				delta_count++;

				//Calculate s matrix
				s[0][0] = s_delta * s_delta * s_delta;
				s[0][1] = s_delta * s_delta;
				s[0][2] = s_delta;
				s[0][3] = 1;

				double q = 0;
				q = findQ(s, H, C);
				interpolate[of][i] = q;
				of++;

			}

			//add last point of each frame to output
			interpolate[of][i] = position[j + 1][i];
			delta_count = 2;
			of++;

		}


	}

	//assume it rest at zero velocity
	for (int i = 0; i < frameLinks; i++) {
		interpolate[oFrames - 1][i] = position[iFrames - 1][i];

	}

	//Output results to robot.ang
	ofstream ofile;
	ofile.open("robot.ang");
	ofile << oFrames;
	ofile << "\n";

	for (int k = 0; k < oFrames; k++) {
		for (int l = 0; l < frameLinks; l++) {
			ofile << interpolate[k][l] << " ";
		}
		ofile << "\n";
	}
	ofile.close();
}

////////////////////////////////////////////////////////Main function//////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {

	/*/Check if exactly 3 Arguments are provided
	if (argc != 2) {
		cout << "Requires robot.key as input" << endl;
		return -1;
	}*/

	//Process the text file and save the points
	parseInputTextFile(argv);

	//Interpolate the points
	positionInterpolation();

	return 0;

}
//////////////////////////////////////////////////////End of Program////////////////////////////////////////////////////////
