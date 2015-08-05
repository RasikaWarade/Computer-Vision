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
#include <math.h>

/////////////////////////////////////////////////////Define Structures////////////////////////////////////////////////////////////////////////////////////

/*
 * Struct: vect for storing x,y,z co-ordinates
 */
struct vect {
	double x;
	double y;
	double z;
};

/*
 *Struct: quaternion to store x,y,z and s co-ordinates for quaternion
 */
struct quaternion {
	double x;
	double y;
	double z;
	double s;
};

//////////////////////////////////////////////////Declare Global variables///////////////////////////////////////////////////////////////////////////////////////

#define MAX_CONFIG  100
using namespace std;
int iFrames, oFrames; //Numer of frames in Input and output frames
std::vector<std::string> words;
double n[1000][3], o[1000][3], a[1000][3], rORp[1000][3]; //Arrays to store n,o and a values for each input frame
double rP[1000][3]; //Array for storing p/r's for output frames
quaternion *quatInput = new quaternion[65530]; //Quaternions for input Frames
quaternion *quatOutput = new quaternion[65530]; //Quaternions for output Frames

////////////////////////////////////////////////////Define Functions///////////////////////////////////////////////////////////////////////////////////

/*
 * Func: to parse the input text file and store the values in the variables
 */
void parseInputTextFile(char* argv[]) {

	bool readFrames = true;

	//text file
	std::string file = "object.key"; //***
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

					//take rest frames
					double v1, v2, v3, v4;
					ss >> v1 >> v2 >> v3 >> v4;
					n[p_cnt][0] = v1;
					o[p_cnt][0] = v2;
					a[p_cnt][0] = v3;
					rORp[p_cnt][0] = v4;

					getline(fin, str);
					std::stringstream ss2;
					ss2 << str;

					ss2 >> v1 >> v2 >> v3 >> v4;
					n[p_cnt][1] = v1;
					o[p_cnt][1] = v2;
					a[p_cnt][1] = v3;
					rORp[p_cnt][1] = v4;
					getline(fin, str);
					std::stringstream ss3;
					ss3 << str;

					ss3 >> v1 >> v2 >> v3 >> v4;
					n[p_cnt][2] = v1;
					o[p_cnt][2] = v2;
					a[p_cnt][2] = v3;
					rORp[p_cnt][2] = v4;
					p_cnt++;
				}
			}
		}

	}
	fin.close();
}

/*
 * Func: calculating Dot product
 */
double dot(vect v1, vect v2) {
	double d;

	d = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	;

	return (d);
}

/*
 * Func: Addition of two quaternions
 */
quaternion Qadd(quaternion q1, quaternion q2) {
	quaternion q;

	q.x = q1.x + q2.x;
	q.y = q1.y + q2.y;
	q.z = q1.z + q2.z;
	q.s = q1.s + q2.s;

	return (q);
}

/*
 *Func: Calculating cross products between two vectors
 */
vect cross(vect v1, vect v2)

{
	vect v;

	v.x = v1.y * v2.z - v1.z * v2.y;
	v.y = v1.z * v2.x - v1.x * v2.z;
	v.z = v1.x * v2.y - v1.y * v2.x;

	return (v);
}

/*
 * Func: Calculating  dot product of Q's
 */

double Qdot(quaternion q1, quaternion q2)

{
	double d;

	d = q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.s * q2.s;

	return (d);
}

/*
 *  Func: find maximum of four numbers
 */
int max_find(double a, double b, double c, double d) {
	int m;

	if (a >= b && a >= c && a >= d)
		m = 1;
	if (b >= a && b >= c && b >= d)
		m = 2;
	if (c >= a && c >= b && c >= d)
		m = 3;
	if (d >= a && d >= b && d >= c)
		m = 4;

	return (m);
}

/*
 * Func: Multiplication of two Q's
 */
quaternion Qmult(quaternion q1, quaternion q2) {
	quaternion q;
	//*** vector cross();
	vect v;
	vect v1, v2;
	v1.x = q1.x;
	v1.y = q1.y;
	v1.z = q1.z;
	v2.x = q1.x;
	v2.y = q2.y;
	v2.z = q2.z;
	v = cross(v1, v2);
	q.x = v.x;
	q.y = v.y;
	q.z = v.z;
	q.x = q.x + q2.s * q1.x + q1.s * q2.x;
	q.y = q.y + q2.s * q1.y + q1.s * q2.y;
	q.z = q.z + q2.s * q1.z + q1.s * q2.z;
	q.s = q1.s * q2.s - dot(v1, v2);

	return (q);
}

/*
 *  Func: Normalized Quaternion
 */
double Qnorm2(quaternion q) {
	double d;					//dot();***
	vect v;
	vect v1;
	v1.x = q.x;
	v1.y = q.y;
	v1.z = q.z;

	d = q.s * q.s + dot(v1, v1);

	return (d);
}

/*
 *  Func: Invert the Quaternion
 */
quaternion Qinv(quaternion q)

{
	quaternion qinv;
	double d;					//,Qnorm2();***

	d = Qnorm2(q);

	qinv.x = -q.x / d;
	qinv.y = -q.y / d;
	qinv.z = -q.z / d;
	qinv.s = q.s / d;

	return (qinv);
}

/*
 * Func: Double of two quaternions
 */
quaternion Double(quaternion q1, quaternion q2) {
	quaternion q;
	double d2;					//***Qdot();

	d2 = 2 * Qdot(q1, q2);

	q.x = d2 * q2.x - q1.x;
	q.y = d2 * q2.y - q1.y;
	q.z = d2 * q2.z - q1.z;
	q.s = d2 * q2.s - q1.s;

	return (q);
}

/*
 * Func: Bisect between the Q's
 */
quaternion Bisect(quaternion q1, quaternion q2) {
	quaternion q;					//***, Qadd();
	double d;					//***, Qnorm2();

	q = Qadd(q1, q2);
	d = sqrt(Qnorm2(q));

	q.x = q.x / d;
	q.y = q.y / d;
	q.z = q.z / d;
	q.s = q.s / d;

	return (q);
}

/*
 * Func: Do spherical interpolation between two Q's
 */
quaternion slerp(quaternion q1, quaternion q2, double u) {

	quaternion qu;
	double sine, cosine, theta, a1, a2;					//*** Qdot();

	cosine = Qdot(q1, q2);
	sine = sqrt(1 - cosine * cosine);
	theta = atan2(sine, cosine);

	if (sine < 0.0001) {
		a1 = 1 - u;
		a2 = u;
	} else {
		a1 = sin((1 - u) * theta) / sine;
		a2 = sin(u * theta) / sine;
	}

	qu.x = a1 * q1.x + a2 * q2.x;
	qu.y = a1 * q1.y + a2 * q2.y;
	qu.z = a1 * q1.z + a2 * q2.z;
	qu.s = a1 * q1.s + a2 * q2.s;

	return (qu);
}

/*
 *Func :  quaternion conversion to rotation matrix conversion
 */
void QtoR(quaternion v) {

	vect n, o, a;
	n.x = 1 - 2 * (v.y * v.y + v.z * v.z);
	n.y = 2 * (v.x * v.y + v.s * v.z);
	n.z = 2 * (v.x * v.z - v.s * v.y);

	o.x = 2 * (v.x * v.y - v.s * v.z);
	o.y = 1 - 2 * (v.x * v.x + v.z * v.z);
	o.z = 2 * (v.y * v.z + v.s * v.x);

	a.x = 2 * (v.x * v.z + v.s * v.y);
	a.y = 2 * (v.y * v.z - v.s * v.x);
	a.z = 1 - 2 * (v.x * v.x + v.y * v.y);

}
/*
 * Func: rotation matrix to quaternion conversion and vice versa
 */

quaternion RtoQ(vect n, vect o, vect a) {

	double dx, dy, dz, ds, s;
	quaternion q;
	vect v;

	dx = 1 + n.x - o.y - a.z;
	dy = 1 - n.x + o.y - a.z;
	dz = 1 - n.x - o.y + a.z;
	ds = 1 + n.x + o.y + a.z;

	switch (max_find(dx, dy, dz, ds)) {
	case 1:
		v.x = sqrt(dx) / 2;
		v.y = (n.y + o.x) / (4 * v.x);
		v.z = (n.z + a.x) / (4 * v.x);
		s = (o.z - a.y) / (4 * v.x);

		break;

	case 2:
		v.y = sqrt(dy) / 2;
		v.x = (n.y + o.x) / (4 * v.y);
		v.z = (o.z + a.y) / (4 * v.y);
		s = (a.x - n.z) / (4 * v.y);

		break;

	case 3:
		v.z = sqrt(dz) / 2;
		v.x = (n.z + a.x) / (4 * v.z);
		v.y = (o.z + a.y) / (4 * v.z);
		s = (n.y - o.x) / (4 * v.z);

		break;

	case 4:
		s = sqrt(ds) / 2;
		v.x = (o.z - a.y) / (4 * s);
		v.y = (a.x - n.z) / (4 * s);
		v.z = (n.y - o.x) / (4 * s);

	}
	q.x = v.x;
	q.y = v.y;
	q.z = v.z;
	q.s = s;
	return q;
}

/*
 * Func: Interpolate
 */
void Qint(int numb_config, quaternion q[]) {
	quaternion a[MAX_CONFIG], b[MAX_CONFIG], qu;
	int i;
	double u, du;

	/* compute intermediate configuration */
	for (i = 1; i < (numb_config - 1); ++i) {

		a[i] = Bisect(Double(q[i - 1], q[i]), q[i + 1]);
		a[i] = slerp(q[i], a[i], 1.0 / 3.0);

		b[i] = Double(a[i], q[i]);

	}

	for (i = 1; i < (numb_config - 2); ++i) {

		for (u = 0; u < 1; u = u + du) {
			qu = slerp(a[i], b[i + 1], u);
			qu = slerp(slerp(slerp(q[i], a[i], u), qu, u),
					slerp(qu, slerp(b[i + 1], q[i + 1], u), u), u);

		}
	}
}

/*
 * Func: find quaternion
 */
double findQ(double s[1][4], double H[4][4], double C[4][1]) {
	double first[4];
	double t = 0;
	///////Matrix Multiplication////////////
	//Find P(u)
	for (int r = 0; r < 4; r++) {	        // 1x4 multiply 4x4
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
 * Func: Perform Orientation Interpolation
 */
void orientationInterpolation() {

	vect n1, o1, a1;
	quaternion q;
	quaternion q1, q2;
	int i = 0;
	for (int k = 0; k < iFrames; k++) {

		n1.x = n[k][0];
		n1.y = n[k][1];
		n1.z = n[k][2];
		o1.x = o[k][0];
		o1.y = o[k][1];
		o1.z = o[k][2];
		a1.x = a[k][0];
		a1.y = a[k][1];
		a1.z = a[k][2];
		q = RtoQ(n1, o1, a1);

		quatInput[i] = q;
		i++;

	}
	//calculate delta
	double delta = (double) ((iFrames - 1)) / (oFrames - 1);
	int noCurves = iFrames - 1;
	int framesMinus = oFrames - iFrames;
	int forEachCurve = (framesMinus - 2) / noCurves;

	int delta_count = 2;
	int of = 2;
	int frames = 0;

	quatOutput[0] = quatInput[0];
	quatOutput[1] = quatInput[0];

	for (int k = 0; k < iFrames - 1; k++) {

		//for each curve
		q1 = quatInput[k];
		q2 = quatInput[k + 1];

		for (int l = 0; l < forEachCurve; l++) {
			double du;
			du = delta * delta_count;
			delta_count++;
			quatOutput[of] = slerp(q1, q2, du);
			of++;
			frames++;
		}
		quatOutput[of] = q2;
		frames++;
		of++;
		delta_count = 2;

	}
	quatOutput[of] = quatInput[iFrames - 1];

	//////////////////////////////////////////////////////////////////////////////////////////
	//for p/r's --input rORp[][3]
	//for each x,y,z calculate final frame px,py,pz

	//Zero frame
	rP[0][0] = rORp[0][0];
	rP[0][1] = rORp[0][1];
	rP[0][2] = rORp[0][2];

	//1st frame
	rP[1][0] = rORp[0][0];
	rP[1][1] = rORp[0][1];
	rP[1][2] = rORp[0][2];

	//Create catmull matrix
	double H[4][4] = { { -0.5, 1.5, -1.5, 0.5 }, // row 0
			{ 1, -2.5, 2, -0.5 }, // row 1
			{ -0.5, 0, 0.5, 0 }, // row 2
			{ 0, 1, 0, 0 } //row 3
	};

	double s[1][4] = { 1, 1, 1, 1 };
	double C1[4][1] = { { 1 }, { 1 }, { 1 }, { 1 } };
	double C2[4][1] = { { 1 }, { 1 }, { 1 }, { 1 } };
	double C3[4][1] = { { 1 }, { 1 }, { 1 }, { 1 } };

	delta_count = 2;
	of = 2;
	frames = 2;

	int cnt = 0;
	double  px3,  py3,  pz3;

	for (int k = 0; k < iFrames - 1; k++) {

		if (k == 0) {

			//x co-ord
			C1[0][0] = rORp[k + 1][0];
			C1[1][0] = rORp[k][0];
			C1[2][0] = rORp[k + 1][0];
			C1[3][0] = rORp[k + 2][0];

			//y co-ord
			C2[0][0] = rORp[k + 1][1];
			C2[1][0] = rORp[k][1];
			C2[2][0] = rORp[k + 1][1];
			C2[3][0] = rORp[k + 2][1];

			//y co-ord
			C3[0][0] = rORp[k + 1][2];
			C3[1][0] = rORp[k][2];
			C3[2][0] = rORp[k + 1][2];
			C3[3][0] = rORp[k + 2][2];

		} else if (k == iFrames - 2) {
			//x co-ord
			C1[0][0] = rORp[k - 1][0];
			C1[1][0] = rORp[k][0];
			C1[2][0] = rORp[k + 1][0];
			C1[3][0] = rORp[k + 1][0];

			//y co-ord
			C2[0][0] = rORp[k - 1][1];
			C2[1][0] = rORp[k][1];
			C2[2][0] = rORp[k + 1][1];
			C2[3][0] = rORp[k + 1][1];

			//y co-ord
			C3[0][0] = rORp[k - 1][2];
			C3[1][0] = rORp[k][2];
			C3[2][0] = rORp[k + 1][2];
			C3[3][0] = rORp[k + 1][2];

		} else if (k == iFrames - 1) {
			//x co-ord
			C1[0][0] = rORp[k - 1][0];
			C1[1][0] = rORp[k][0];
			C1[2][0] = rORp[k][0];
			C1[3][0] = rORp[k][0];

			//y co-ord
			C2[0][0] = rORp[k - 1][1];
			C2[1][0] = rORp[k][1];
			C2[2][0] = rORp[k][1];
			C2[3][0] = rORp[k][1];

			//y co-ord
			C3[0][0] = rORp[k - 1][2];
			C3[1][0] = rORp[k][2];
			C3[2][0] = rORp[k][2];
			C3[3][0] = rORp[k][2];
		} else {
			//x co-ord
			C1[0][0] = rORp[k - 1][0];
			C1[1][0] = rORp[k][0];
			C1[2][0] = rORp[k + 1][0];
			C1[3][0] = rORp[k + 2][0];

			//y co-ord
			C2[0][0] = rORp[k - 1][1];
			C2[1][0] = rORp[k][1];
			C2[2][0] = rORp[k + 1][1];
			C2[3][0] = rORp[k + 2][1];

			//y co-ord
			C3[0][0] = rORp[k - 1][2];
			C3[1][0] = rORp[k][2];
			C3[2][0] = rORp[k + 1][2];
			C3[3][0] = rORp[k + 2][2];
		}
		for (int l = 0; l < forEachCurve; l++) {
			double s_delta = 0;
			s_delta = delta * delta_count;
			delta_count++;

			//Calculate s matrix
			s[0][0] = s_delta * s_delta * s_delta;
			s[0][1] = s_delta * s_delta;
			s[0][2] = s_delta;
			s[0][3] = 1;

			px3 = findQ(s, H, C1);
			py3 = findQ(s, H, C2);
			pz3 = findQ(s, H, C3);

			rP[of][0] = px3;
			rP[of][1] = py3;
			rP[of][2] = pz3;

			of++;
			frames++;
		}

		rP[of][0] = rORp[k + 1][0];
		rP[of][1] = rORp[k + 1][1];
		rP[of][2] = rORp[k + 1][2];

		of++;
		frames++;

		delta_count = 2;

	}

	//last
	rP[oFrames - 1][0] = rORp[iFrames - 1][0];
	rP[oFrames - 1][1] = rORp[iFrames - 1][1];
	rP[oFrames - 1][2] = rORp[iFrames - 1][2];

	//////////////////////////////////////////////////////////////////////////////////////////

	//Output results to object.traj
	ofstream ofile;
	ofile.open("object.traj");
	ofile << oFrames;
	ofile << "\n";

	for (int i = 0; i < oFrames; i++) {

		//QtoR
		vect n, o, a;
		n.x = 1
				- 2
						* (quatOutput[i].y * quatOutput[i].y
								+ quatOutput[i].z * quatOutput[i].z);
		n.y = 2
				* (quatOutput[i].x * quatOutput[i].y
						+ quatOutput[i].s * quatOutput[i].z);
		n.z = 2
				* (quatOutput[i].x * quatOutput[i].z
						- quatOutput[i].s * quatOutput[i].y);

		o.x = 2
				* (quatOutput[i].x * quatOutput[i].y
						- quatOutput[i].s * quatOutput[i].z);
		o.y = 1
				- 2
						* (quatOutput[i].x * quatOutput[i].x
								+ quatOutput[i].z * quatOutput[i].z);
		o.z = 2
				* (quatOutput[i].y * quatOutput[i].z
						+ quatOutput[i].s * quatOutput[i].x);

		a.x = 2
				* (quatOutput[i].x * quatOutput[i].z
						+ quatOutput[i].s * quatOutput[i].y);
		a.y = 2
				* (quatOutput[i].y * quatOutput[i].z
						- quatOutput[i].s * quatOutput[i].x);
		a.z = 1
				- 2
						* (quatOutput[i].x * quatOutput[i].x
								+ quatOutput[i].y * quatOutput[i].y);

		ofile << n.x << " " << o.x << " " << a.x << " " << rP[i][0] << "\n";
		ofile << n.y << " " << o.y << " " << a.y << " " << rP[i][1] << "\n";
		ofile << n.z << " " << o.z << " " << a.z << " " << rP[i][2] << "\n";

	}

	ofile.close();
}

///////////////////////////////////////////////////Main Function///////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
	/*
	//Check if exactly 3 Arguments are provided
	if (argc != 2) {
		cout << "Requires object.key as input" << endl;
		return -1;
	}*/

	//Process the text file and save the points
	parseInputTextFile(argv);

	//Interpolate the points
	orientationInterpolation();

	return 0;

}
/////////////////////////////////////////////////End of Program/////////////////////////////////////////////////////////
