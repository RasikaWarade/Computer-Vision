#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <cstring>
#include <string>
#include<math.h>
#include<algorithm>
#include <locale> 
#include <string>
#include <algorithm>
#include <cctype> 
using namespace std;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

template<typename T>
inline T& max(T& x, T& y) {
	return x > y ? x : y;
}

template<typename U>
inline U& min(U& x, U& y) {
	return x < y ? x : y;
}

//Data structure for Ray,Sphere,Camera,Light,Vertices,Shade
struct Camera {
	string camera_name;
	float p_x, p_y, p_z, v_x, v_y, v_z;
	float vup_x, vup_y, vup_z;
	float near, far;
};

struct Ray {
	string scene_cast;
	int w, h, recurr;
};

struct Light {
	float x, y, z, w;
	double Ipx, Ipy, Ipz;
};

struct Sphere {
	string sphere_name;
	float X, Y, Z;
	int R, G, B;
	float radius;
	float v, c, b, dsquare, s;
	string usemtl;
};

struct Vect {
	float x, y, z;
};

struct Shade {
	string usemtl;
	float alpa;
	double ka_r, ka_g, ka_b, kd_r, kd_g, kd_b, ks_r, ks_g, ks_b;
	double n1, Tr, Kr, Krf;
};

struct Face {
	string usemtl;
	string name;
	float a, b, c;

};

struct ToUpper {
	char operator()(char c) const {
		return std::toupper(c);
	}
};

struct RGB {
	double r, g, b;
};
Vect *vect = new Vect[65530];

int c_no, s_no, r_no, l_no, m_no, v_no, f_no;

int img_w = 0, img_h = 0;
int c_pixel[6000][6000][3];
int depth_pixel[6000][6000][3];
Camera *cam = new Camera[65530];
Ray *scene = new Ray[65530];
Light *lights = new Light[65530];
Shade *material = new Shade[65530];
Vect *vertices = new Vect[65530];
Sphere *sph = new Sphere[65530];
Face *face = new Face[65530];
Vect *vert = new Vect[65530];
float SMIN, TMIN;
int R, G, B;
void fileParsing() {
	c_no = l_no = s_no = r_no = m_no = v_no = f_no = 0;

	// create a file-reading object
	ifstream fin;
	std::vector<std::string> parsed_token;
	std::string str;
	std::string color;
	std::string filename;
	//object file
	std::string str8;
	//str8=argv[1];
	str8 = "";
	//std::cout << "Splitting: " << str8 << '\n';
	unsigned found = str8.find_last_of("/\\");
	//std::cout << " path: " << str8.substr(0,found) << '\n';
	//std::cout << " file: " << str8.substr(found+1) << '\n';
	std::string path = str8.substr(0, found);
	std::string filep;
	filep = str8.substr(found + 1);
	//cout<<filep;

	//fin.open(filep.c_str());
	 //fin.open(filep.c_str(),std::ifstream::in);//**for now opening from the same folder
	fin.open("model.obj");
	if (!fin.good())
		cout << "Could not open file\n";
	else {
		// read each line of the file
		while (!fin.eof()) {
			getline(fin, str);

			stringstream ss;
			ss << str;
			string s;
			int flag = 0;
			ss >> s;
			std::transform(s.begin(), s.end(), s.begin(), ToUpper());

			// cout<<s<<endl;
			if (!s.compare("MTLLIB")) {  //***
				ss >> filename;

				//cout<<color<<endl;

			} else if (!s.compare("USEMTL")) {
				ss >> color;
				flag = 0;

			} else if (!s.compare("S")) {
				if (flag == 0) {
					sph[s_no].usemtl = color;
					flag = 1;
				}
				ss >> sph[s_no].sphere_name >> sph[s_no].X >> sph[s_no].Y
						>> sph[s_no].Z >> sph[s_no].radius;
				s_no++;
			} else if (!s.compare("V")) {

				ss >> vertices[v_no].x >> vertices[v_no].y >> vertices[v_no].z;
				v_no++;
			} else if (!s.compare("G")) {
				ss >> face[f_no].name;
				face[1].usemtl = "NA";

			} else if (!s.compare("F")) {

				int number_of_words = 0;
				int word;
				int f[50];
				int i = 0;

				if (flag == 0) {
					face[f_no].usemtl = color;
					flag = 1;
				}
				while (ss >> word) {
					f[i] = word;
					number_of_words++;
					i++;
				}

				//need to generalise for as many number of vertices in faces***
				if (number_of_words == 4) {
					face[f_no].a = f[0];
					face[f_no].b = f[1];
					face[f_no].c = f[2];
					f_no++;
					face[f_no].a = f[0];
					face[f_no].b = f[2];
					face[f_no].c = f[3];
					f_no++;

				} else if (number_of_words == 3) {
					face[f_no].a = f[0];
					face[f_no].b = f[1];
					face[f_no].c = f[2];
					f_no++;
				}

			}

		}
	}
	fin.close();
	//object file ends

	//Commands file
	//str8=argv[2];
	str8 = "";
	//std::cout << "Splitting: " << str8 << '\n';
	unsigned found1 = str8.find_last_of("/\\");
	//std::cout << " path: " << str8.substr(0,found1) << '\n';
	//std::cout << " file: " << str8.substr(found1+1) << '\n';
	path = str8.substr(0, found1);

	filep = str8.substr(found1 + 1);
	//cout<<filep;
	// fin.open(filep.c_str(),std::ifstream::in);//**for now opening from the same folder
	fin.open("commands");
	//fin.open(filep.c_str());
	if (!fin.good())
		cout << "Could not open file\n";
	else {
		// read each line of the file
		while (!fin.eof()) {
			getline(fin, str);

			stringstream ss;
			ss << str;
			string s;
			ss >> s;
			std::transform(s.begin(), s.end(), s.begin(), ToUpper());

			if (!s.compare("C")) {
				ss >> cam[c_no].camera_name >> cam[c_no].p_x >> cam[c_no].p_y
						>> cam[c_no].p_z >> cam[c_no].v_x >> cam[c_no].v_y
						>> cam[c_no].v_z >> cam[c_no].vup_x >> cam[c_no].vup_y
						>> cam[c_no].vup_z >> cam[c_no].near >> cam[c_no].far;

				c_no++;
			} else if (!s.compare("L")) {
				ss >> lights[l_no].x >> lights[l_no].y >> lights[l_no].z
						>> lights[l_no].w >> lights[l_no].Ipx
						>> lights[l_no].Ipy >> lights[l_no].Ipz;

				l_no++;
			} else if (!s.compare("R")) {
				ss >> scene[r_no].scene_cast >> scene[r_no].w >> scene[r_no].h
						>> scene[r_no].recurr;
				r_no++;
			}

		}
	}
	fin.close();
	//commands file done
	 filename=path+"/"+filename;
		//fin.open(filename.c_str());
		//fin.open(filename.c_str(),std::ifstream::in);

	//fin.open(filename.c_str());
	fin.open("model.mtl");
	if (!fin.good())
		cout << "Could not open file\n";
	else {
		// read each line of the file

		while (!fin.eof()) {
			getline(fin, str);

			stringstream ss;
			ss << str;
			string s;
			string ill;
			int flag = 0;
			ss >> s;
			std::transform(s.begin(), s.end(), s.begin(), ToUpper());

			if (!s.compare("NEWMTL")) {
				ss >> material[m_no].usemtl;
				flag = 1;
			} else if (!s.compare("KA")) {
				ss >> material[m_no].ka_r >> material[m_no].ka_g
						>> material[m_no].ka_b;

			} else if (!s.compare("KD")) {
				ss >> material[m_no].kd_r >> material[m_no].kd_g
						>> material[m_no].kd_b;
			} else if (!s.compare("KS")) {
				ss >> material[m_no].ks_r >> material[m_no].ks_g
						>> material[m_no].ks_b;
			} else if (!s.compare("NS")) {
				ss >> material[m_no].alpa;
				//m_no++;
			} else if (!s.compare("N1")) {
				ss >> material[m_no].n1;
				//m_no++;
			} else if (!s.compare("TR")) {
				ss >> material[m_no].Tr;
				//m_no++;
			} else if (!s.compare("KR")) {
				ss >> material[m_no].Kr;
				//m_no++;
			} else if (!s.compare("KRF")) {
				ss >> material[m_no].Krf;
				m_no++;
			}
		}

	}
	fin.close();
	//mtlfile ends*/
}

void saveToPPM(int e, int f) {
	//e-for each ..f-for each camera
	//Color image
	ofstream color_ppm;
	std::ostringstream color;
	color << scene[e].scene_cast << "_" << cam[f].camera_name << "_color.ppm";
	std::string var = color.str();
	color_ppm.open(var.c_str());

	color_ppm << "P3 " << img_w << " " << img_h << " " << "255 \n" << endl;

	for (int l = 0; l < scene[e].w; l++) {
		for (int k = 0; k < scene[e].h; k++) {
			color_ppm << c_pixel[k][img_h - 1 - l][0] << " "
					<< c_pixel[k][img_h - 1 - l][1] << " "
					<< c_pixel[k][img_h - 1 - l][2] << "\n" << endl;
		}
	}
	color_ppm.close();

	//Depth image***depth for every object
	ofstream depth_ppm;
	std::ostringstream dpt;
	dpt << scene[e].scene_cast << "_" << cam[f].camera_name << "_depth.ppm";
	std::string var1 = dpt.str();

	depth_ppm.open(var1.c_str());

	depth_ppm << "P3 " << img_w << " " << img_h << " " << "255 \n" << endl;

	for (int l = 0; l < scene[e].w; l++) {
		for (int k = 0; k < scene[e].h; k++) {
			depth_ppm << depth_pixel[k][img_h - 1 - l][0] << " "
					<< depth_pixel[k][img_h - 1 - l][1] << " "
					<< depth_pixel[k][img_h - 1 - l][2] << "\n" << endl;
		}
	}
	depth_ppm.close();

}

int findShadow(float actual_x, float actual_y, float actual_z, float Lx,
		float Ly, float Lz, float m4) {

	if (s_no != 0)
		for (int t = 0; t < s_no; t++)	//for each object
				{
			//if(index!=t)
			//{
			float n_c_x, n_c_y, n_c_z;
			n_c_x = sph[t].X - actual_x;
			n_c_y = sph[t].Y - actual_y;
			n_c_z = sph[t].Z - actual_z;

			//Ray vector
			float v, b, dsquare, c;
			float s_s;
			float n_v_x, n_v_y, n_v_z;
			//Ray direction
			n_v_x = Lx * n_c_x;
			n_v_y = Ly * n_c_y;
			n_v_z = Lz * n_c_z;
			v = n_v_x + n_v_y + n_v_z;

			c = sqrt(
					((sph[t].X - actual_x) * (sph[t].X - actual_x))
							+ ((sph[t].Y - actual_y) * (sph[t].Y - actual_y))
							+ ((sph[t].Z - actual_z) * (sph[t].Z - actual_z)));

			b = sqrt((c * c) - (v * v));

			dsquare = (sph[t].radius * sph[t].radius) - (b * b);

			if (dsquare > 0) {
				s_s = v - (sqrt(dsquare));
				if (s_s < m4 && s_s > 0.00001) {		//*****
				//cout<<" pixel in shadow"<<endl;
				//isInShadow=true;
					return 1;
				}
			}
			//}
		}
	//if(nearest!=1000)//**bet 0 and 1 then shadow
	//return 1;
	return 0;
}

Vect computeReflectionRay(Vect V, Vect N) {
	Vect R;
	//Generate Reflection Ray
	//R=V-2(V.N)N
	float dot;

	dot = V.x * N.x + V.y * N.y + V.z * N.z;
	if (dot > 0) {
		R.x = (2 * dot * N.x) - V.x;
		R.y = (2 * dot * N.y) - V.y;
		R.z = (2 * dot * N.z) - V.z;

		//cout<<R.x<<" "<<R.y<<" "<<R.z<<endl;

	}
	return R;
}

RGB Local(float actual_x, float actual_y, float actual_z, float N_x, float N_y,
		float N_z, float V_x, float V_y, float V_z, int colo) {
	RGB local;
	//for each light
	//float L_mag;
	float L_x, L_y, L_z;
	int shadow = 0;
	double I_r, I_g, I_b;
	double Ia = 20;
	//float actual_x,actual_y,actual_z;
	//float N_x,N_y,N_z,L_x,L_y,L_z;
	//float N_mag;
	float L_mag;
	float diffuse, total_diffuse;
	float zero = 0;
	//double range=255;
	//float V_x,V_y,V_z;
	//float m2,m3;
	float R_x, R_y, R_z;
	float V_dot_R;
	float specular;
	//cout<<l_no<<endl;
	//Ambient reflection
	I_r = Ia * material[colo].ka_r;
	I_g = Ia * material[colo].ka_g;
	I_b = Ia * material[colo].ka_b;

	for (int ll = 0; ll < l_no; ll++) {

		L_mag =
				sqrt(
						((lights[ll].x - actual_x) * (lights[ll].x - actual_x))
								+ ((lights[ll].y - actual_y)
										* (lights[ll].y - actual_y))
								+ ((lights[ll].z - actual_z)
										* (lights[ll].z - actual_z)));
		L_x = (lights[ll].x - actual_x) / L_mag;
		L_y = (lights[ll].y - actual_y) / L_mag;
		L_z = (lights[ll].z - actual_z) / L_mag;

		//Dot product of N and L
		//check dot product value and apply conditions

		shadow = findShadow(actual_x, actual_y, actual_z, L_x, L_y, L_z, L_mag);
		//shadow=0;
		//if background color
		//intersect light with each sphere
		//PixelP0  Light(Lx,Ly,Lz)

		if (shadow) {
			//dnt add light
			//I_r=0;
			//I_g=0;
			//I_b=0;
		} else {
			diffuse = (N_x * L_x) + (N_y * L_y) + (N_z * L_z);

			//if(diffuse<0)
			//diffuse=diffuse*(-1);

			diffuse = max(diffuse, zero);
			//cout<<"diffue "<<diffuse<<endl;

			I_r = I_r + (lights[ll].Ipx * material[colo].kd_r * diffuse);
			I_g = I_g + (lights[ll].Ipy * material[colo].kd_g * diffuse);
			I_b = I_b + (lights[ll].Ipz * material[colo].kd_b * diffuse);

			//Calculate Specular highlight

			//Calculate R=2(L.N)N-L
			R_x = (2 * diffuse * N_x) - L_x;
			R_y = (2 * diffuse * N_y) - L_y;
			R_z = (2 * diffuse * N_z) - L_z;

			//Dot product of V and R
			V_dot_R = (V_x * R_x) + (V_y * R_y) + (V_z * R_z);
			V_dot_R = max(V_dot_R, zero);

			specular = pow(V_dot_R, material[colo].alpa);

			I_r = I_r + (lights[ll].Ipx * material[colo].ks_r * specular);
			I_g = I_g + (lights[ll].Ipy * material[colo].ks_g * specular);
			I_b = I_b + (lights[ll].Ipz * material[colo].ks_b * specular);

		}
	}
	local.r = I_r;
	local.g = I_g;
	local.b = I_b;
	return local;
}

int checkSphereIntersection(Vect org, Vect dir) {
	int index = -1;
	float smin = 10000;
	if (s_no != 0)		//check if sphere objects exists in scene
		for (int i = 0; i < s_no; i++) {
			float rsquare;
			rsquare = sph[i].radius * sph[i].radius;
			float c_x, c_y, c_z;
			c_x = sph[i].X - org.x;
			c_y = sph[i].Y - org.y;
			c_z = sph[i].Z - org.z;
			float c_mag = sqrt((c_x * c_x) + (c_y * c_y) + (c_z * c_z));
			float v_x, v_y, v_z;
			float v;					// c, b, dsquare, s;
			v_x = dir.x * c_x;
			v_y = dir.y * c_y;
			v_z = dir.z * c_z;
			v = v_x + v_y + v_z;
			//float v_mag = sqrt(v*v);
			float dsquare = rsquare - ((c_mag * c_mag) - (v * v));
			//cout<<"dsquarre "<<dsquare<<" cmag "<<c_mag <<" v_mag "<<v<<" rsquare "<<rsquare<<endl;
			if (dsquare >= 0) {
				float d = sqrt(dsquare);
				float t1 = v - d;
				if (t1 > 0) {
					if (t1 < smin) {
						//cout<<"smin "<<smin<<endl;
						smin = t1;
						index = i;
						sph[i].s = smin;

					}

				}
			}

		}
	//cout<<"Shphere index "<<index<<" smin "<<smin<<" SMIN "<<SMIN<<endl;

	if (index != -1) {
		SMIN = smin;
		return index;

	} else {
		SMIN = smin;
		return -1;
	}

}
int checkPolygonIntersection(Vect org, Vect dir) {
	float beta, gamma, tt;
	int index = -1;
	float smin = 10000;
	if (f_no != 0)					//check if sphere objects exists in scene
		for (int jj = 0; jj < f_no; jj++) {
			int v0, v1, v2;
			v0 = face[jj].a - 1;
			v1 = face[jj].b - 1;
			v2 = face[jj].c - 1;
			//Get the three point of triangle
			vert[0].x = vertices[v0].x;
			vert[0].y = vertices[v0].y;
			vert[0].z = vertices[v0].z;

			vert[1].x = vertices[v1].x;
			vert[1].y = vertices[v1].y;
			vert[1].z = vertices[v1].z;

			vert[2].x = vertices[v2].x;
			vert[2].y = vertices[v2].y;
			vert[2].z = vertices[v2].z;

			//3 euations in 3 unknowns
			//L=>PRP
			//U=>x1-PRP/|x1-PRP|
			//

			float aa, bb, cc, dd, ee, ff, gg, hh, ii, pp, qq, rr;
			//(B-A)===>A(a,b,c)
			aa = vert[3].x = vert[1].x - vert[0].x;
			dd = vert[3].y = vert[1].y - vert[0].y;
			gg = vert[3].z = vert[1].z - vert[0].z;
			//(C-A)==>B(d,e,f)
			bb = vert[4].x = vert[2].x - vert[0].x;
			ee = vert[4].y = vert[2].y - vert[0].y;
			hh = vert[4].z = vert[2].z - vert[0].z;
			//C==>-U(g,h,i)
			float Ux, Uy, Uz;
			float mag;
			//mag=sqrt(((x1-org.x)*(x1-org.x))+((y1-org.y)*(y1-org.y))+((z1-org.z)*(z1-org.z)));
			cc = Ux = -dir.x;							//(x1-cam[f].p_x)/mag;
			ff = Uy = -dir.y;							//(y1-cam[f].p_y)/mag;
			ii = Uz = -dir.z;							//(z1-cam[f].p_z)/mag;

			//D==>L-A(p,q,r)
			float Dx, Dy, Dz;
			pp = Dx = org.x - vert[0].x;
			qq = Dy = org.y - vert[0].y;
			rr = Dz = org.z - vert[0].z;

			//find beta, gamma and t

			float dtr = aa * ee * ii + bb * ff * gg + cc * dd * hh
					- (cc * ee * gg + aa * ff * hh + dd * bb * ii);
			float aaa = ee * ii - ff * hh;
			float bbb = cc * hh - bb * ii;
			float ccc = bb * ff - cc * ee;
			float ddd = ff * gg - dd * ii;
			float eee = aa * ii - cc * gg;
			float fff = cc * dd - aa * ff;
			float ggg = dd * hh - ee * gg;
			float hhh = bb * gg - aa * hh;
			float iii = aa * ee - bb * dd;
			if (dtr != 0) {
				beta = (pp * aaa + qq * bbb + rr * ccc) / dtr;
				gamma = (pp * ddd + qq * eee + rr * fff) / dtr;
				tt = (pp * ggg + qq * hhh + rr * iii) / dtr;
			} else {
				beta = 0;
				gamma = 0;
				tt = 0;

			}

			if (beta >= 0 && gamma >= 0 && (beta + gamma) <= 1) {
				//float d = sqrt(dsquare);
				//float t1 = v - d;
				if (tt > 0) {
					if (tt < smin) {
						//cout<<"smin "<<smin<<endl;
						smin = tt;
						index = jj;
						//sph[i].s = smin;

					}

				}
			}

		}
	//cout<<"Poly index "<<index<<" smin "<<smin<<" TMIN "<<TMIN<<endl;

	if (index != -1) {
		TMIN = smin;
		return index;

	} else {
		TMIN = smin;
		return -1;
	}
}

RGB GetColor(Vect org, Vect dir, int depth) {
	float range = 255;
	RGB intensity, local, spec;
	RGB phong;
	RGB irefl;
	RGB irefr;
	int inter_obj = -1;
	int index = -1;
	int index_p = -1;

	if (depth == 0) {

		intensity.r = 0;
		intensity.b = 0;
		intensity.g = 0;
		return intensity;
	}

	else {
		//check if the ray intersects any object
		//cout<<"depth "<<depth<<endl;
		SMIN = 10000;
		index = checkSphereIntersection(org, dir);

		//cout<<inter_obj<<endl;
		TMIN = 10000;
		index_p = checkPolygonIntersection(org, dir);
		//Check polygon intervention here

		if (index == -1 && index_p == -1) {
			intensity.r = 0;
			intensity.b = 0;
			intensity.g = 0;
			return intensity;
		} else {
			//First-assign any color to the sphere
			//compare SMIN AND TMIN
			//in the minimum one assing intersection_obj
			//Calculate Shadow and phong for polygon and sphere
			//local=Phong(org,dir,inter_obj);
			//intensity=local;

			if (SMIN > TMIN) {
				float actual_x, actual_y, actual_z;
				float poly_r = 0, poly_g = 0, poly_b = 0;
				float beta, gamma, tt;
				int index4 = 10000;

				float tmin = 10000;
				if (f_no != 0)			//Check for each polygonal face
					for (int jj = 0; jj < f_no; jj++) {

						int v0, v1, v2;
						v0 = face[jj].a - 1;
						v1 = face[jj].b - 1;
						v2 = face[jj].c - 1;
						//Get the three point of triangle
						vert[0].x = vertices[v0].x;
						vert[0].y = vertices[v0].y;
						vert[0].z = vertices[v0].z;

						vert[1].x = vertices[v1].x;
						vert[1].y = vertices[v1].y;
						vert[1].z = vertices[v1].z;

						vert[2].x = vertices[v2].x;
						vert[2].y = vertices[v2].y;
						vert[2].z = vertices[v2].z;

						//3 euations in 3 unknowns
						//L=>PRP
						//U=>x1-PRP/|x1-PRP|
						//

						float aa, bb, cc, dd, ee, ff, gg, hh, ii, pp, qq, rr;
						//(B-A)===>A(a,b,c)
						aa = vert[3].x = vert[1].x - vert[0].x;
						dd = vert[3].y = vert[1].y - vert[0].y;
						gg = vert[3].z = vert[1].z - vert[0].z;
						//(C-A)==>B(d,e,f)
						bb = vert[4].x = vert[2].x - vert[0].x;
						ee = vert[4].y = vert[2].y - vert[0].y;
						hh = vert[4].z = vert[2].z - vert[0].z;
						//C==>-U(g,h,i)
						float Ux, Uy, Uz;
						float mag;
						//mag=sqrt(((x1-cam[f].p_x)*(x1-cam[f].p_x))+((y1-cam[f].p_y)*(y1-cam[f].p_y))+((z1-cam[f].p_z)*(z1-cam[f].p_z)));
						cc = Ux = -dir.x;				//(x1-cam[f].p_x)/mag;
						ff = Uy = -dir.y;				//(y1-cam[f].p_y)/mag;
						ii = Uz = -dir.z;

						//D==>L-A(p,q,r)
						float Dx, Dy, Dz;
						pp = Dx = org.x - vert[0].x;
						qq = Dy = org.y - vert[0].y;
						rr = Dz = org.z - vert[0].z;

						//find beta, gamma and t

						float dtr = aa * ee * ii + bb * ff * gg + cc * dd * hh
								- (cc * ee * gg + aa * ff * hh + dd * bb * ii);
						float aaa = ee * ii - ff * hh;
						float bbb = cc * hh - bb * ii;
						float ccc = bb * ff - cc * ee;
						float ddd = ff * gg - dd * ii;
						float eee = aa * ii - cc * gg;
						float fff = cc * dd - aa * ff;
						float ggg = dd * hh - ee * gg;
						float hhh = bb * gg - aa * hh;
						float iii = aa * ee - bb * dd;
						if (dtr != 0) {
							beta = (pp * aaa + qq * bbb + rr * ccc) / dtr;
							gamma = (pp * ddd + qq * eee + rr * fff) / dtr;
							tt = (pp * ggg + qq * hhh + rr * iii) / dtr;
						} else {
							beta = 0;
							gamma = 0;
							tt = 0;

						}
						if (tt > 0) {
							if (beta >= 0 && gamma >= 0
									&& (beta + gamma) <= 1) {
								if (tt < tmin) {
									tmin = tt;
									index_p = jj;
								}
							}

						}
					}
						if (index_p !=-1) {
							inter_obj = index_p;

							vect[5].x = vect[3].x;
							vect[5].y = vect[3].y;
							vect[5].z = vect[3].z;

							vect[6].x = vect[4].x;
							vect[6].y = vect[4].y;
							vect[6].z = vect[4].z;

							double I_r, I_g, I_b;
							double Ia = 20;
							int count;
							int colo = 0;
							if (!face[1].usemtl.compare("NA")) {
								count = 0;
							} else
								count = index_p;

							while ((face[count].usemtl).compare(
									material[colo].usemtl))	//***change it for each face
								colo++;
							//cout<<face[0].usemtl<<endl;
							//cout<<"Material "<<material[colo].usemtl<<" "<<material[colo].ks_r<<" "<<material[colo].ks_g<<" "<<material[colo].ks_b<<" "<<material[colo].kd_r<<" "<<material[colo].kd_g<<" "<<material[colo].kd_b<<" "<<material[colo].ka_r<<" "<<material[colo].ka_g<<" "<<material[colo].ka_b<<endl;
							//Ambient reflection
							I_r = Ia * material[colo].ka_r;
							I_g = Ia * material[colo].ka_g;
							I_b = Ia * material[colo].ka_b;

							//Diffuse Reflection

							//Ip*Kd*(N.L)

							float N_x, N_y, N_z, L_x, L_y, L_z;
							float N_mag;
							float L_mag;
							float diffuse, total_diffuse;
							float zero = 0;
							double range = 255;
							float V_x, V_y, V_z;
							float m2, m3;
							float R_x, R_y, R_z;
							float V_dot_R;
							float specular;
							float diff_r = 0, diff_b = 0, diff_g = 0,
									spec_r = 0, spec_g = 0, spec_b = 0;
							int shadow = 0;
							//find co-ordinates of	the point (x,y,z) on the sphere
							//PRP to (x1,y1,z1) plus s gives (x,y,z)

							actual_x = org.x + (tmin * dir.x);
							actual_y = org.y + (tmin * dir.y);
							actual_z = org.z + (tmin * dir.z);
							//cout<<tmin<<endl;
							//cout<<actual_x<<" "<<actual_y<<" "<<actual_z<<endl;
							//For polygon N=E1xE2/|E1xE2|

							/*
							 N_mag=sqrt(((actual_x-sph[index].X)*(actual_x-sph[index].X))+((actual_y-sph[index].Y)*(actual_y-sph[index].Y))+((actual_z-sph[index].Z)*(actual_z-sph[index].Z)));
							 N_x=(actual_x-sph[index].X)/N_mag;
							 N_y=(actual_y-sph[index].Y)/N_mag;
							 N_z=(actual_z-sph[index].Z)/N_mag;//Surface Normal
							 * */

							////

							//Vert[3]=E1 Vert[4]=E2
							N_x = vert[3].y * vert[4].z - vert[4].y * vert[3].z;
							N_y = vert[4].x * vert[3].z - vert[3].x * vert[4].z;
							N_z = vert[3].x * vert[4].y - vert[3].y * vert[4].x;
							N_mag = sqrt(
									(N_x * N_x) + (N_y * N_y) + (N_z * N_z));
							N_x = N_x / N_mag;
							N_y = N_y / N_mag;
							N_z = N_z / N_mag;

							//calculate V=C-Q/|C-Q|

							m2 = sqrt(
									((org.x - actual_x) * (org.x - actual_x))
											+ ((org.y - actual_y)
													* (org.y - actual_y))
											+ ((org.z - actual_z)
													* (org.z - actual_z)));
							V_x = (org.x - actual_x) / m2;
							V_y = (org.y - actual_y) / m2;
							V_z = (org.z - actual_z) / m2;
							//cout<<" V_x"<<V_x<<" V_y "<<V_y<<" V_z "<< V_z<<endl;

							//Specular
							Vect R;
							Vect PRP;
							PRP.x = actual_x;
							PRP.y = actual_y;
							PRP.z = actual_z;
							//R = computeReflectionRay(refl_V, refl_N);
							//R=V-2(V.N)N

							float dot;

							dot = V_x * N_x + V_y * N_y + V_z * N_z;
							//if(dot>0){
							R.x = (2 * dot * N_x) - V_x;
							R.y = (2 * dot * N_y) - V_y;
							R.z = (2 * dot * N_z) - V_z;
							irefl = GetColor(PRP, R, depth - 1);

							//For each light
							for (int ll = 0; ll < l_no; ll++) {
								float N_L;
								L_mag = sqrt(
										((lights[ll].x - actual_x)
												* (lights[ll].x - actual_x))
												+ ((lights[ll].y - actual_y)
														* (lights[ll].y
																- actual_y))
												+ ((lights[ll].z - actual_z)
														* (lights[ll].z
																- actual_z)));
								L_x = (lights[ll].x - actual_x) / L_mag;
								L_y = (lights[ll].y - actual_y) / L_mag;
								L_z = (lights[ll].z - actual_z) / L_mag;

								//Dot product of N and L
								//check dot product value and apply conditions
								//Dot product of N and L
								//check dot product value and apply conditions
								N_L = L_x * N_x + L_y * N_y + L_z * N_z;
								if (N_L >= 0) {
									//if the pixel is in shadow

									shadow = findShadow(actual_x, actual_y,
											actual_z, L_x, L_y, L_z, L_mag);
									//shadow=0;
									//if background color
									//intersect light with each sphere
									//PixelP0  Light(Lx,Ly,Lz)

									if (shadow) {
										//dnt add light
									} else {
										diffuse = (N_x * L_x) + (N_y * L_y)
												+ (N_z * L_z);

										//if(diffuse<0)
										//diffuse=diffuse*(-1);

										diffuse = max(diffuse, zero);
										//cout<<"diffue "<<diffuse<<endl;

										I_r = I_r
												+ (lights[ll].Ipx
														* material[colo].kd_r
														* diffuse);
										I_g = I_g
												+ (lights[ll].Ipy
														* material[colo].kd_g
														* diffuse);
										I_b = I_b
												+ (lights[ll].Ipz
														* material[colo].kd_b
														* diffuse);

										//Calculate Specular highlight

										//Calculate R=2(L.N)N-L
										R_x = (2 * diffuse * N_x) - L_x;
										R_y = (2 * diffuse * N_y) - L_y;
										R_z = (2 * diffuse * N_z) - L_z;

										//Dot product of V and R
										V_dot_R = (V_x * R_x) + (V_y * R_y)
												+ (V_z * R_z);
										V_dot_R = max(V_dot_R, zero);

										specular = pow(V_dot_R,
												material[colo].alpa);

										I_r = I_r
												+ (lights[ll].Ipx
														* material[colo].ks_r
														* specular);
										I_g = I_g
												+ (lights[ll].Ipy
														* material[colo].ks_g
														* specular);
										I_b = I_b
												+ (lights[ll].Ipz
														* material[colo].ks_b
														* specular);

									}

									I_r=min(I_r,range);
									I_g=min(I_g,range);
									I_b=min(I_b,range);

									poly_r = I_r;
									poly_g = I_g;
									poly_b = I_b;
								}
							}

							///2
						}



				local.r = poly_r + 1 * irefl.r;
				local.g = poly_g + 1 * irefl.g;
				local.b = poly_b + 1 * irefl.b;

				intensity = local;
				//intensity.r=0;
				//intensity.b=255;
				//intensity.g=0;
			} else if (SMIN < TMIN) {
				//local=Phong(org,dir,index,depth);
				int colo = 0;
				while ((sph[index].usemtl).compare(material[colo].usemtl))
					colo++;

				//cout<<"Material "<<material[colo].usemtl<<" "<<material[colo].ks_r<<" "<<material[colo].ks_g<<" "<<material[colo].ks_b<<" "<<material[colo].kd_r<<" "<<material[colo].kd_g<<" "<<material[colo].kd_b<<" "<<material[colo].ka_r<<" "<<material[colo].ka_g<<" "<<material[colo].ka_b<<endl;

				//Diffuse Reflection

				//Ip*Kd*(N.L)

				float actual_x, actual_y, actual_z;
				float N_x, N_y, N_z;// L_x, L_y, L_z;
				float N_mag;
				float L_mag;
				float diffuse, total_diffuse;
				float zero = 0;
				double range=255;
				float V_x, V_y, V_z;
				float m2;//, m3;
				float R_x, R_y, R_z;
				float V_dot_R;
				float specular;

				//int shadow = 0;
				//float diff_r = 0, diff_b = 0, diff_g = 0, spec_r = 0,
					//	spec_g = 0, spec_b = 0;

				//THIS IS ACTUAL HIT
				//       org						dir
				actual_x = org.x + (sph[index].s * dir.x);//(x1-cam[f].p_x))/m;
				actual_y = org.y + (sph[index].s * dir.y);//(y1-cam[f].p_y))/m;
				actual_z = org.z + (sph[index].s * dir.z);//(z1-cam[f].p_z))/m;

				//For Sphere N=(S-Pc)/|S-Pc|
				//L=P-S/|P-S|
				N_mag = sqrt(
						((actual_x - sph[index].X) * (actual_x - sph[index].X))
								+ ((actual_y - sph[index].Y)
										* (actual_y - sph[index].Y))
								+ ((actual_z - sph[index].Z)
										* (actual_z - sph[index].Z)));
				N_x = (actual_x - sph[index].X) / N_mag;
				N_y = (actual_y - sph[index].Y) / N_mag;
				N_z = (actual_z - sph[index].Z) / N_mag;		//Surface Normal

				//calculate V=C-Q/|C-Q|

				m2 = sqrt(
						((org.x - actual_x) * (org.x - actual_x))
								+ ((org.y - actual_y) * (org.y - actual_y))
								+ ((org.z - actual_z) * (org.z - actual_z)));
				V_x = (org.x - actual_x) / m2;
				V_y = (org.y - actual_y) / m2;
				V_z = (org.z - actual_z) / m2;
				//cout<<" V_x"<<V_x<<" V_y "<<V_y<<" V_z "<< V_z<<endl;

				//for each light
				//Check shadow and calculate phong
				phong = Local(actual_x, actual_y, actual_z, N_x, N_y, N_z, V_x,
						V_y, V_z, colo);
				//Calculate reflected ray
				Vect refl_V;
				// =V_x;				//W=signed opposite of the ray being traced
				//V_y =-V_y;
				//V_z =-V_z;

				////Vect refl_N;
				//refl_N.x = N_x;
				//refl_N.y = N_y;
				//refl_N.z = N_z;
				//Calculate the reflected component
				RGB irefr;
				//cout<<material[colo].Tr<<endl;
			{
				Vect R;
				Vect PRP;
				PRP.x = actual_x;
				PRP.y = actual_y;
				PRP.z = actual_z;
				//R = computeReflectionRay(refl_V, refl_N);
				//R=V-2(V.N)N

				float dot;

				dot = V_x * N_x + V_y * N_y + V_z * N_z;
				//if(dot>0){
				R.x = (2 * dot * N_x) - V_x;
				R.y = (2 * dot * N_y) - V_y;
				R.z = (2 * dot * N_z) - V_z;
				irefl = GetColor(PRP, R, depth - 1);


				}
			if(material[colo].Tr<1){

					//Refract
					//cout<<"inside refraction"<<endl;

					Vect T;
					Vect refr;

					//refr=computeRefractionRay(Vect refr,Vect N);

					float currentInd=material[colo].n1;
					float airInd=1.0;
					float u=airInd/currentInd;

					float a=0,b=0;
					float wn;
					//cout<<"u"<<u<<endl;
					wn = V_x * N_x + V_y * N_y + V_z * N_z;
					//cout<<"wn "<<wn<<endl;
					a=-u;

					b=sqrt((((u*u)*(wn*wn))-(u*u)+1));
					//cout<<"a "<<a<<" b "<<b<<endl;

					if(b>=0){
						b=a*wn+sqrt((((u*u)*(wn*wn))-(u*u)+1));
					T.x=(a*V_x)+(b*N_x);
					T.y=(a*V_y)+(b*N_y);
					T.z=(a*V_z)+(b*N_z);
					//cout<<T.x<<" "<<T.y<<" "<<" "<<T.z<<endl;

					//Normalize T
					float t_mag=0;
					t_mag=sqrt((T.x*T.x)+(T.y*T.y)+(T.z*T.z));
					T.x=T.x/t_mag;
					T.y=T.y/t_mag;
					T.z=T.z/t_mag;
					//cout<<T.x<<" "<<T.y<<" "<<" "<<T.z<<endl;
					//cout<<t_mag<<endl;
					//Got T
					//Now Get Refr ray
					//float v_dist=0;
					float vv_x,vv_y,vv_z,vv;

					vv_x=(sph[index].X-actual_x)*T.x;
					vv_y=(sph[index].Y-actual_y)*T.y;
					vv_z=(sph[index].Z-actual_z)*T.z;
					vv=vv_x+vv_y+vv_z;

					if(vv>0)
					{
					//cout<<"vv" <<vv<<" radius "<<sph[index].radius<<endl;
					//find the next hit point inside sphere
					//cout<<"ACTUAL "<<actual_x<<" "<<actual_y<<" "<<actual_z<<endl;
					Vect next;
					next.x=actual_x+(2*vv*T.x);
					next.y=actual_y+(2*vv*T.y);
					next.z=actual_z+(2*vv*T.z);
					//cout<<" NEXT "<<next.x<<" "<<next.y<<" "<<next.z<<endl;
					//Find the normal again--check the sign

					float n_m;
					Vect Norm;
					n_m=sqrt(((sph[index].X-next.x)*(sph[index].X-next.x))+((sph[index].Y-next.y)*(sph[index].Y-next.y))+((sph[index].Z-next.z)*(sph[index].Z-next.z)));
									Norm.x = (sph[index].X-next.x) / n_m;
									Norm.y = (sph[index].Y-next.y) / n_m;
									Norm.z = (sph[index].Z-next.z) / n_m;		//Surface Normal
					Vect V_2;
					V_2.x=-T.x;
					V_2.y=-T.y;
					V_2.z=-T.z;

					/////

					float currentIn = material[colo].n1;
						float nextInd = 1.0;
						float u2 = currentIn / nextInd;

						float a2 = 0, b2 = 0;
						float wn2;
						wn2 = V_2.x * Norm.x + V_2.y * Norm.y + V_2.z * Norm.z;
						a2 = -u2;

						b2 =  sqrt((((u2 * u2) * (wn2 * wn2)) - (u2 * u2) + 1));
						if (b2 >= 0) {
							b2 = a2 * wn2 + sqrt((((u2 * u2) * (wn2 * wn2)) - (u2 * u2) + 1));
							refr.x = (a2 * V_2.x) + (b2 * Norm.x);
							refr.y = (a2 * V_2.y) + (b2 * Norm.y);
							refr.z = (a2 * V_2.z) + (b2 * Norm.z);
							//Normalize T
							float refr_mag = 0;
							refr_mag = sqrt(
									(refr.x * refr.x) + (refr.y * refr.y) + (refr.z * refr.z));
							refr.x = refr.x / refr_mag;
							refr.y = refr.y / refr_mag;
							refr.z = refr.z / refr_mag;

							irefr = GetColor(next, refr, depth - 1);
							//cout<<"irefr "<<irefr.r<<" "<<irefr.g<<" "<<irefr.b<<endl;


						}
					}
					}
				}

				local.r =(material[colo].Tr*(phong.r + material[colo].Kr * irefl.r))+((1-material[colo].Tr)*irefr.r);
				local.g =(material[colo].Tr*(phong.g + material[colo].Kr * irefl.g))+((1-material[colo].Tr)*irefr.g);
				local.b =(material[colo].Tr*(phong.b + material[colo].Kr * irefl.b))+((1-material[colo].Tr)*irefr.b);
				intensity = local;

				//intensity.r=255;
				//intensity.b=255;
				//intensity.g=0;
			}

		}
	}
	return intensity;
}

void trace(int e, int f) {

	float d = cam[f].near;
	//Centre of image plane-->
	vect[0].x = cam[f].p_x - (d * cam[f].v_x);
	vect[0].y = cam[f].p_y - (d * cam[f].v_y);
	vect[0].z = cam[f].p_z - (d * cam[f].v_z);

	//N-->Vect[1]	
	float n_mag;
	vect[1].x = cam[f].v_x;
	vect[1].y = cam[f].v_y;
	vect[1].z = cam[f].v_z;
	n_mag = sqrt(
			((cam[f].v_x) * (cam[f].v_x)) + ((cam[f].v_y) * (cam[f].v_y))
					+ ((cam[f].v_z) * (cam[f].v_z)));

	if (n_mag != 0) {
		vect[1].x = cam[f].v_x / n_mag;
		vect[1].y = cam[f].v_y / n_mag;
		vect[1].z = cam[f].v_z / n_mag;
	}
	//U-->Vect[2]
	vect[2].x = cam[f].vup_y * vect[1].z - vect[1].y * cam[f].vup_z;
	vect[2].y = vect[1].x * cam[f].vup_z - cam[f].vup_x * vect[1].z;
	vect[2].z = cam[f].vup_x * vect[1].y - cam[f].vup_y * vect[1].x;

	float u_mag;
	u_mag = sqrt(
			((vect[2].x) * (vect[2].x)) + ((vect[2].y) * (vect[2].y))
					+ ((vect[2].z) * (vect[2].z)));

	if (u_mag != 0) {
		vect[2].x = vect[2].x / u_mag;
		vect[2].y = vect[2].y / u_mag;
		vect[2].z = vect[2].z / u_mag;
	}
	//V-->Vect[3]
	vect[3].x = vect[1].y * vect[2].z - vect[2].y * vect[1].z;
	vect[3].y = vect[2].x * vect[1].z - vect[1].x * vect[2].z;
	vect[3].z = vect[1].x * vect[2].y - vect[1].y * vect[2].x;

	//cout<<"| "<<vect[2].x<<" "<<vect[2].y<<" "<<vect[2].z<<" 0 |"<<endl;
	//cout<<"| "<<vect[3].x<<" "<<vect[3].y<<" "<<vect[3].z<<" 0 |"<<endl;
	//cout<<"| "<<vect[1].x<<" "<<vect[1].y<<" "<<vect[1].z<<" 0 |"<<endl;
	//cout<<"| 0"<<" 0"<<" 0"<<" 1 |"<<endl;

	float u_img, v_img, m;
	//scene[e].w=50;
	//scene[e].h=50;
	///FOR IMAGE
	for (int x = 0; x < scene[e].w; x++)
		for (int y = 0; y < scene[e].h; y++) {
			int w = scene[e].w - 1;
			int h = scene[e].h - 1;

			//Compute image plane co-ordinated bounded in (-1,-1) to (1,1)
			u_img = float(x) * 2 / w - 1;
			v_img = float(y) * 2 / h - 1;

			float x1, y1, z1;
			x1 = cam[f].p_x + (u_img * vect[2].x) + (v_img * vect[3].x)
					- (d * vect[1].x);
			y1 = cam[f].p_y + (u_img * vect[2].y) + (v_img * vect[3].y)
					- (d * vect[1].y);
			z1 = cam[f].p_z + (u_img * vect[2].z) + (v_img * vect[3].z)
					- (d * vect[1].z);

			m = sqrt(
					((x1 - cam[f].p_x) * (x1 - cam[f].p_x))
							+ ((y1 - cam[f].p_y) * (y1 - cam[f].p_y))
							+ ((z1 - cam[f].p_z) * (z1 - cam[f].p_z)));

			//initialize IMAGE to default
			c_pixel[x][y][0] = 0;
			c_pixel[x][y][1] = 0;
			c_pixel[x][y][2] = 0;

			depth_pixel[x][y][0] = depth_pixel[x][y][1] = depth_pixel[x][y][2] =
					0;

			//CONTRUCT RAY THROUGH EACH PIXEL
			Vect dir;
			dir.x = ((x1 - cam[f].p_x)) / m;
			dir.y = ((y1 - cam[f].p_y)) / m;
			dir.z = ((z1 - cam[f].p_z)) / m;

			Vect org;
			org.x = cam[f].p_x;
			org.y = cam[f].p_y;
			org.z = cam[f].p_z;

			//int maxRecurr=scene[e].recurr;
			int maxRecurr = scene[e].recurr;
			//maxRecurr=1;
			//GET THE COLOR FOR THE PIXEL
			RGB color;
			color = GetColor(org, dir, maxRecurr);
			double range = 255;
			double zero = 0;
			color.r = max(color.r, zero);
			color.g = max(color.g, zero);
			color.b = max(color.b, zero);

			color.r = min(color.r, range);
			color.g = min(color.g, range);
			color.b = min(color.b, range);

			//FINAL GLOBAL ILLUMINATION FOR THE PIXEL
			c_pixel[x][y][0] = (int) color.r;
			c_pixel[x][y][1] = (int) color.g;
			c_pixel[x][y][2] = (int) color.b;

		}

	//Output the color_values and depth_values to save in ppm
	img_w = scene[e].w;
	img_h = scene[e].h;
	saveToPPM(e, f);
}

int main(int argc,char *argv[]) {
	//fileparsing(argv);
	fileParsing();
	/*for(int i=0;i<m_no;i++)
	 {

	 cout<<material[i].usemtl<<" " <<material[i].ka_r<<" "<<material[i].ka_g<<" "<<material[i].ka_b<<" "<<material[i].kd_r<<" "<<material[i].kd_g<<" "<<material[i].kd_b<<" "<<material[i].ks_r<<" "<<material[i].ks_g<<" "<<material[i].ks_b<<" "<<material[i].alpa<<material[i].n1<<" "<<material[i].Tr<<" "<<material[i].Kr<<" "<<material[i].Krf<<endl;
	 }
	 /*
	 //sphere
	 for(int i=0;i<s_no;i++)
	 {

	 cout<<sph[i].usemtl<<" "<<sph[i].sphere_name<<" " <<sph[i].X<<" "<<sph[i].Y<<" "<<sph[i].Z<<" "<<sph[i].radius<<endl;
	 }
	 //vertices
	 for(int i=0;i<v_no;i++)
	 {

	 cout<<vertices[i].x<<vertices[i].y<<vertices[i].z<<endl;
	 }
	 for(int i=0;i<f_no;i++)
	 {

	 cout<<face[i].a<<face[i].b<<face[i].c<<endl;
	 cout<<face[i].usemtl<<endl;
	 }*/
	for (int a = 0; a < r_no; a++) {

		//for each scene
		for (int b = 0; b < c_no; b++) {
			//for each camera
			trace(a, b);

		}
	}
	return 0;
}
