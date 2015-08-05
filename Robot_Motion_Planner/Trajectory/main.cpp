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
#include<math.h>

/////////////////////////////////////////////////////Declare Global variable////////////////////////////////////////////////////////////////////////////////////

typedef struct {
    float x;
    float y;

} VECT;

#define MAXDOF  10    /* maximum number of degrees of freedom for robot */
#define LAMBDA2  0.021999    /* damping factor squared */

float theta[MAXDOF]; /* joint angle positions */
float thetaOver[MAXDOF]; /* joint angle velocities */

float rr[MAXDOF];
float lambda[MAXDOF];

VECT efjacob[MAXDOF]; /* Jacobian for end effector */
VECT ediff[MAXDOF];

using namespace std;
int noOfJoints, noOfPositions;
float norm = 0.0;

float linkLength[1000][2], position[1000][2];

////////////////////////////////////////////////////Define Functions///////////////////////////////////////////////////////////////////////////////////
/*
 *Calculate end effector Jacobian
 */
void efjacobian(int no) {
    int i, j;
    float psi; /* absolute joint angle */
    VECT aj[MAXDOF]; /* absolute angle jacobian */

    /* calculate terms of Jacobian in absolute angles */
    psi = 0;
    for (i = 0; i < noOfJoints; i++) {
        psi = thetaOver[i];
        aj[i].x = -linkLength[i][0] * sin(psi);
        aj[i].y = linkLength[i][0] * cos(psi);
    }

    /* Jacobian in relative angles is a sum of absolute terms */
    //cout<<"JACOBIAN IN TERMS OF theta"<<endl;
    for (i = 0; i < noOfJoints; i++) {
        efjacob[i].x = 0;
        efjacob[i].y = 0;
        for (j = i; j < noOfJoints; j++) {
            efjacob[i].x = efjacob[i].x + aj[j].x;
            efjacob[i].y = efjacob[i].y + aj[j].y;
        }
        //cout << efjacob[i].x << " " << efjacob[i].y << endl;
    }
}

void solve(int no) /* solve the inverse kinematics using damped least squares */
{
    int i;
    float determ, temp;
    VECT efpseudo[MAXDOF];

    float jjT[2][2];

    /* calculate the Jacobian times the Jacobian transpose (jjT) */
    jjT[0][0] = 0;
    jjT[1][0] = 0;
    jjT[0][1] = 0;
    jjT[1][1] = 0;

    for (i = 0; i < noOfJoints; i++) {
        jjT[0][0] = jjT[0][0] + efjacob[i].x * efjacob[i].x;
        jjT[1][0] = jjT[1][0] + efjacob[i].y * efjacob[i].x;

        jjT[0][1] = jjT[0][1] + efjacob[i].x * efjacob[i].y;
        jjT[1][1] = jjT[1][1] + efjacob[i].y * efjacob[i].y;

    }

    /* add damping factor to the diagonal */
    jjT[0][0] = jjT[0][0] +lambda[no];//+ LAMBDA2;
    jjT[1][1] = jjT[1][1] +lambda[no];//+ LAMBDA2;

    /* calculate the inverse of jjT */
    determ = jjT[0][0] * jjT[1][1] - jjT[0][1] * jjT[1][0];
    temp = jjT[0][0];
    jjT[0][0] = jjT[1][1] / determ;
    jjT[0][1] = -jjT[0][1] / determ;
    jjT[1][0] = -jjT[1][0] / determ;
    jjT[1][1] = temp / determ;

    //cout << "Jt inverse..." << endl;
    //cout << jjT[0][0] << " " << jjT[0][1] << endl;
    //cout << jjT[1][0] << " " << jjT[1][1] << endl;

    //cout<<"MULTIPLY BY JT inverse"<<endl;
    /* multiply by jT for damped least squares inverse, i.e. jT(jjT)^-1 */
    for (i = 0; i < noOfJoints; i++) {
        efpseudo[i].x = 0;
        efpseudo[i].y = 0;
        efpseudo[i].x = efpseudo[i].x + efjacob[i].x * jjT[0][0]
                + efjacob[i].y * jjT[1][0];
        efpseudo[i].y = efpseudo[i].y + efjacob[i].x * jjT[0][1]
                + efjacob[i].y * jjT[1][1];
        //cout << efpseudo[i].x << " " << efpseudo[i].y << endl;
    }

    //9 terms
    //cout<<"Edifferece"<<" "<<ediff[no].x<< " "<<ediff[no].y<<endl;
    //matrix multiply
    for (i = 0; i < noOfJoints; i++) {
        theta[i] = 0;

        theta[i] = efpseudo[i].x * ediff[no].x + efpseudo[i].y * ediff[no].y;

    }
}
/*
 * Func: Parse the input text file and store the values in the variables
 */

void parseInputTextFile(char* argv[]) {
    //REad arm file
    //text file
    std::string file = "arm";

    //create a file-reading object
    ifstream fin;
    fin.open(file.c_str(), std::ifstream::in);

    //Check if file is open
    if (!fin.good())
        cout << "Could not open arm file\n";
    else {
        std::string str;

        //Read file and save the points
        int p_cnt = 0;
        while (!fin.eof()) {
            getline(fin, str);
            std::stringstream ss;
            ss << str;
            if (str.length() != 0) {
                if (p_cnt == 0) {
                    //read the noofjoints and norm
                    ss >> noOfJoints >> norm;
                    p_cnt++;
                } else {
                    //read the linklengths and starting postion
                    ss >> linkLength[p_cnt - 1][0] >> linkLength[p_cnt - 1][1];
                    p_cnt++;

                }

            }
        }

    }
    fin.close();

    //********************************************************//
    //REad trajectory file
    //text file
    std::string file2 = "trajectory";

    //create a file-reading object
    ifstream fin2;
    fin2.open(file2.c_str(), std::ifstream::in);

    //Check if file is open
    if (!fin2.good())
        cout << "Could not open trajectory file\n";
    else {
        std::string str2;

        //Read file and save the points
        int p_cnt2 = 0;
        while (!fin2.eof()) {
            getline(fin2, str2);
            std::stringstream ss2;
            ss2 << str2;
            if (str2.length() != 0) {
                if (p_cnt2 == 0) {
                    //read the noofjoints and norm
                    ss2 >> noOfPositions;
                    p_cnt2++;
                } else {
                    //read the linklengths and starting postion
                    ss2 >> position[p_cnt2 - 1][0] >> position[p_cnt2 - 1][1];
                    p_cnt2++;

                }

            }
        }
    }
    fin2.close();

}

////////////////////////////////////////////////////////Main function//////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {

    //Process the text file and save the points
    parseInputTextFile(argv);

    /*
     cout<<"NoofJoints:"<<noOfJoints<<endl;
     cout<<"NoofPositions:"<<noOfPositions<<endl;
     cout<<"Norm:"<<norm<<endl;

     for(int i=0;i<noOfJoints;i++){
     cout<<linkLength[i][0]<<" "<<linkLength[i][1]<<endl;
     }

     for(int i=0;i<noOfPositions;i++){
     cout<<position[i][0]<<" "<<position[i][1]<<endl;
     }
     */

    //cout << "Calculate end effector velocity..." << endl;
    //ediff[0].x=0;
    //ediff[0].y=0;
    for (int i = 0; i < noOfPositions - 1; i++) {
        ediff[i].x = position[i + 1][0] - position[i][0];
        ediff[i].y = position[i + 1][1] - position[i][1];
        //cout << ediff[i].x << " " << ediff[i].y << endl;
        rr[i]=(ediff[i].x *ediff[i].x +ediff[i].y*ediff[i].y);


        lambda[i]=rr[i]*rr[i];
        lambda[i]=lambda[i]/0.2;
    }

    //Output results to robot.ang
    ofstream ofile;
    ofile.open("angles");

    //###/* convert to radians */

    for (int i = 0; i < noOfPositions; i++) {
        if (i == 0) {
            for (int j = 0; j < noOfJoints; j++) {

                thetaOver[j] = linkLength[j][1];
            }


            for (int j = 0; j < noOfJoints; j++) {
                cout << " " << thetaOver[j];
                ofile<<" "<<thetaOver[j];
            }
            cout << endl;
            ofile<<"\n";

            efjacobian(i);
            //
            solve(i);
            //


        } else {

            for (int j = 0; j < noOfJoints; j++) {

                thetaOver[j] = thetaOver[j] + theta[j];

            }

            efjacobian(i);
            //
            solve(i);


            //cout<<"Theta values.."<<endl;
            for (int j = 0; j < noOfJoints; j++) {
                cout << " " << thetaOver[j];
                ofile<<" "<<thetaOver[j];
            }
            cout << endl;
            ofile<<"\n";
        }
    }

    return 0;

}
//////////////////////////////////////////////////////End of Program////////////////////////////////////////////////////////
