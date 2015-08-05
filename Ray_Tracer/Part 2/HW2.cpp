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
inline T& max(T& x, T& y)
{
return x>y ? x:y;
}

template<typename U>
inline U& min(U& x, U& y)
{
return x<y ? x:y;
}

//Data structure for Ray,Sphere,Camera,Light,Vertices,Shade
struct Camera{
string camera_name;
float p_x,p_y,p_z,v_x,v_y,v_z;
float vup_x,vup_y,vup_z;
float near,far;
};

struct Ray{
string scene_cast;
int w,h,recurr;
};

struct Light{
	float x,y,z,w;
	double Ipx,Ipy,Ipz;
};

struct Sphere{
string sphere_name;
float X,Y,Z;
int R,G,B;
float radius;
float v,c,b,dsquare,s;
string usemtl;
};


struct Vect{
	float x,y,z;
};

struct Shade{
	string usemtl;
	float alpa;
	double ka_r,ka_g,ka_b,kd_r,kd_g,kd_b,ks_r,ks_g,ks_b;
};

struct Face{
	string usemtl;
	string name;
	float a,b,c;

};

 struct ToUpper
   {
     char operator() (char c) const  { return std::toupper(c); }
   };
   
Vect *vect=new Vect[65530];


int c_no,s_no,r_no,l_no,m_no,v_no,f_no;
int img_w=0,img_h=0;
int c_pixel[6000][6000][3];
int depth_pixel[6000][6000][3];
Camera *cam=new Camera[65530];
Ray *scene=new Ray[65530];
Light *lights=new Light[65530];
Shade *material=new Shade[65530];
Vect *vertices=new Vect[65530];
Sphere *sph=new Sphere[65530];
Face *face=new Face[65530];
Vect *vert=new Vect[65530];

int R,G,B;
void fileParsing(char *argv[])
{
	c_no=l_no=s_no=r_no=m_no=v_no=f_no=0;
	
  // create a file-reading object
  ifstream fin;
  std::vector<std::string> parsed_token;
  std::string str;
  std:: string color;
  std::string filename;
  //object file
  std::string str8;
  str8=argv[1];
  //std::cout << "Splitting: " << str8 << '\n';
  unsigned found = str8.find_last_of("/\\");
  //std::cout << " path: " << str8.substr(0,found) << '\n';
  //std::cout << " file: " << str8.substr(found+1) << '\n';
  std::string path= str8.substr(0,found);
  std::string filep;
  filep=str8.substr(found+1);
  //cout<<filep;		
  
  fin.open(filep.c_str());
  if (!fin.good()) 
   cout<<"Could not open file\n";
  else{
  // read each line of the file
  while (!fin.eof())
  {
	 getline(fin,str);

	 stringstream ss;
	 ss<<str;
	 string s;
	 int flag=0;
	 ss>>s;
	 std::transform (s.begin(), s.end(), s.begin(), ToUpper());
	  
	 // cout<<s<<endl;
	  if(!s.compare("MTLLIB")){//***
		ss>>filename;
		
		//cout<<color<<endl;
		
	 } else
	 if(!s.compare("USEMTL")){
		ss>>color;
		flag=0;
			
	 }
	 else if(!s.compare("S")){
	 if(flag==0){
			sph[s_no].usemtl=color;
			flag=1;
	 }
	 ss>>sph[s_no].sphere_name>>sph[s_no].X>>sph[s_no].Y>>sph[s_no].Z>>sph[s_no].radius;
	 s_no++;
	 }
	 else if(!s.compare("V")){
	 
	 ss>>vertices[v_no].x>>vertices[v_no].y>>vertices[v_no].z;
	 v_no++;
	 }
	 else if(!s.compare("G")){
		 ss>>face[f_no].name;
		 face[1].usemtl="NA";
		
	 }
	 else if(!s.compare("F")){
	
	 int number_of_words = 0;
	 int word;
	 int f[50];
	 int i=0;
	 
	 if(flag==0){
			face[f_no].usemtl=color;
			flag=1;
	 }
	 while (ss >> word)
	 {
		 f[i]=word;
		number_of_words++;
		i++;
	 }
	
	 //need to generalise for as many number of vertices in faces***
	 if(number_of_words==4)
	 {
		 face[f_no].a=f[0];
		 face[f_no].b=f[1];
		 face[f_no].c=f[2];
		 f_no++;
		 face[f_no].a=f[0];
		 face[f_no].b=f[2];
		 face[f_no].c=f[3];
		 f_no++;


	 }
	 else if(number_of_words==3){
		  face[f_no].a=f[0];
		 face[f_no].b=f[1];
		 face[f_no].c=f[2];
		 f_no++;
	 }

	 }
	 
	
  }
  }
  fin.close();
  //object file ends

 //Commands file
 
  str8=argv[2];
  //std::cout << "Splitting: " << str8 << '\n';
  unsigned found1 = str8.find_last_of("/\\");
  //std::cout << " path: " << str8.substr(0,found1) << '\n';
  //std::cout << " file: " << str8.substr(found1+1) << '\n';
  path= str8.substr(0,found1);

  filep=str8.substr(found1+1);
  //cout<<filep;		
 // fin.open(argv[2]);
 fin.open(filep.c_str());
  if (!fin.good()) 
   cout<<"Could not open file\n";
  else{
  // read each line of the file
  while (!fin.eof())
  {
	 getline(fin,str);

	 stringstream ss;
	 ss<<str;
	 string s;
	 ss>>s;
	 std::transform (s.begin(), s.end(), s.begin(), ToUpper());
	
	 if(!s.compare("C")){
		 ss>>cam[c_no].camera_name>>cam[c_no].p_x>>cam[c_no].p_y>>cam[c_no].p_z>>cam[c_no].v_x>>cam[c_no].v_y>>cam[c_no].v_z>>cam[c_no].vup_x>>cam[c_no].vup_y>>cam[c_no].vup_z>>cam[c_no].near>>cam[c_no].far;

		 c_no++;
	 }
	 else if (!s.compare("L")){
		 ss>>lights[l_no].x>>lights[l_no].y>>lights[l_no].z>>lights[l_no].w>>lights[l_no].Ipx>>lights[l_no].Ipy>>lights[l_no].Ipz;	

		  l_no++;
	 }
	 else if (!s.compare("R")){
		 ss>>scene[r_no].scene_cast>>scene[r_no].w>>scene[r_no].h>>scene[r_no].recurr;
			r_no++;
	 }
	
  }
  }
  fin.close();
   //commands file done

  
  fin.open(filename.c_str());
  if (!fin.good()) 
   cout<<"Could not open file\n";
  else{
	   // read each line of the file
 
  while (!fin.eof())
  {
	 getline(fin,str);

	 stringstream ss;
	 ss<<str;
	 string s;
	 string ill;
	 int flag=0;
	 ss>>s;
	 std::transform (s.begin(), s.end(), s.begin(), ToUpper());
	
	if(!s.compare("NEWMTL")){
		ss>>material[m_no].usemtl;
		flag=1;		
	}
	 else
	 if(!s.compare("KA")){
		 ss>>material[m_no].ka_r>>material[m_no].ka_g>>material[m_no].ka_b;
		
	 }
	 else
		 if(!s.compare("KD")){
			 ss>>material[m_no].kd_r>>material[m_no].kd_g>>material[m_no].kd_b;
	 }
	 else
		 if(!s.compare("KS")){
			 ss>>material[m_no].ks_r>>material[m_no].ks_g>>material[m_no].ks_b;
	 }
	 else
		 if(!s.compare("NS")){
			 ss>>material[m_no].alpa;
			 m_no++;
		 }
	 }
	
   
  }
  fin.close();
  //mtlfile ends*/
}

void saveToPPM(int e,int f){
	//e-for each ..f-for each camera
	//Color image
	ofstream color_ppm;	
	std::ostringstream color;
	color <<scene[e].scene_cast<<"_"<<cam[f].camera_name<<"_color.ppm";
	std::string var = color.str();
	color_ppm.open(var.c_str());
	
	color_ppm<<"P3 "<< img_w<<" "<<img_h<<" "<<"255 \n"<<endl;
	
	for(int l=0;l<scene[e].w;l++)
	{
		for(int k=0;k<scene[e].h;k++)
		{
			color_ppm<<c_pixel[k][img_h-1-l][0]<<" "<<c_pixel[k][img_h-1-l][1]<<" "<<c_pixel[k][img_h-1-l][2]<<"\n"<<endl;
		}
	}
	color_ppm.close();
	
	//Depth image***depth for every object
	ofstream depth_ppm;
	std::ostringstream dpt;
	dpt <<scene[e].scene_cast<<"_"<<cam[f].camera_name<<"_depth.ppm";
	std::string var1 = dpt.str();
		
	depth_ppm.open(var1.c_str());
	
	depth_ppm<<"P3 "<< img_w<<" "<<img_h<<" "<<"255 \n"<<endl;
	
	for(int l=0;l<scene[e].w;l++)
	{
		for(int k=0;k<scene[e].h;k++)
		{
			depth_ppm<<depth_pixel[k][img_h-1-l][0]<<" "<<depth_pixel[k][img_h-1-l][1]<<" "<<depth_pixel[k][img_h-1-l][2]<<"\n"<<endl;
		}
	}
	depth_ppm.close();
	
}

int findShadow(float actual_x,float actual_y,float actual_z,float Lx,float Ly,float Lz,float m4,float index){
	
	float nearest=1000;
	if(s_no!=0)
	for(int t=0;t<s_no;t++)//for each object
	{
		//if(index!=t)
		//{
		float n_c_x,n_c_y,n_c_z;
		n_c_x=sph[t].X-actual_x;
		n_c_y=sph[t].Y-actual_y;
		n_c_z=sph[t].Z-actual_z;

		//Ray vector
		float n_v_x,n_v_y,n_v_z;
		//Ray direction
		n_v_x=Lx*n_c_x;
		n_v_y=Ly*n_c_y;
		n_v_z=Lz*n_c_z;
		sph[t].v=n_v_x+n_v_y+n_v_z;	
						
		sph[t].c=sqrt(((sph[t].X-actual_x)*(sph[t].X-actual_x))+((sph[t].Y-actual_y)*(sph[t].Y-actual_y))+((sph[t].Z-actual_z)*(sph[t].Z-actual_z)));
			
		sph[t].b=sqrt((sph[t].c*sph[t].c)-(sph[t].v*sph[t].v));

		sph[t].dsquare=(sph[t].radius*sph[t].radius)-(sph[t].b*sph[t].b);
		
		if(sph[t].dsquare>0)
				{	
					sph[t].s=sph[t].v-(sqrt(sph[t].dsquare));
					if(sph[t].s<m4 && sph[t].s>0.00001){//*****
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

void trace(int e,int f)
{
	//cout<<"e"<<e<<"F"<<f<<endl;
	//cout<<"Started Rendering image.."<<endl;
	//Part1: Position and PRP
	//concentrating on camera placement
	//c camera_name prp_x prp_y prp_z vpn_x vpn_y vpn_z vup_x vup_y vup_z near far
	//c cam02 2.5 2.5 0 2.5 2.5 0 0 1 0 1 5
	
	float d=cam[f].near;	
	//Centre of image plane-->
	vect[0].x=cam[f].p_x-(d*cam[f].v_x);
	vect[0].y=cam[f].p_y-(d*cam[f].v_y);
	vect[0].z=cam[f].p_z-(d*cam[f].v_z);
	
	//cout<<"Centre of image plane : c_x  "<<vect[0].x<<"  c_y  "<<vect[0].y<<"  c_z  "<<vect[0].z<<endl;
	
	//N-->Vect[1]	
	float n_mag;
	vect[1].x=cam[f].v_x;
	vect[1].y=cam[f].v_y;
	vect[1].z=cam[f].v_z;
	n_mag=sqrt(((cam[f].v_x)*(cam[f].v_x))+((cam[f].v_y)*(cam[f].v_y))+((cam[f].v_z)*(cam[f].v_z)));
	//cout<<"n_mag "<<n_mag<<endl;
	if(n_mag!=0)
	{
	vect[1].x=cam[f].v_x/n_mag;
	vect[1].y=cam[f].v_y/n_mag;
	vect[1].z=cam[f].v_z/n_mag;
	//cout<<" n_x  "<<vect[1].x<<"  n_y  "<<vect[1].y<<"  n_z  "<<vect[1].z<<endl;
	}
	
	//cout<<"cam2[0].vup_x "<<cam2[0].vup_x<<" cam2[0].vup_y "<<cam2[0].vup_y<<" cam2[0].vup_z "<<cam2[0].vup_z<<endl;
	//cout<<" vect[1].x "<<vect[1].x<<" vect[1].y "<<vect[1].y<<" vect[1].z "<<vect[1].z<<endl;

	//U-->Vect[2]
	vect[2].x = cam[f].vup_y*vect[1].z - vect[1].y*cam[f].vup_z;	
	vect[2].y =  vect[1].x*cam[f].vup_z-cam[f].vup_x*vect[1].z;
	vect[2].z  = cam[f].vup_x*vect[1].y - cam[f].vup_y*vect[1].x; 
	
	float u_mag;
	u_mag=sqrt(((vect[2].x)*(vect[2].x))+((vect[2].y)*(vect[2].y))+((vect[2].z)*(vect[2].z)));
	//cout<<"u_mag "<<u_mag<<endl;
	if(u_mag!=0){
	vect[2].x=vect[2].x/u_mag;
	vect[2].y=vect[2].y/u_mag;
	vect[2].z=vect[2].z/u_mag;
	//cout<<" U_x  "<<vect[2].x<<"  U_y  "<<vect[2].y<<"  U_z  "<<vect[2].z<<endl;
	}
	
	//V-->Vect[3]
	
	vect[3].x = vect[1].y*vect[2].z - vect[2].y*vect[1].z;
	vect[3].y = vect[2].x*vect[1].z - vect[1].x*vect[2].z;
	vect[3].z  = vect[1].x*vect[2].y - vect[1].y*vect[2].x; 
	
	
	//cout<<" V_x  "<<vect[3].x<<"  V_y  "<<vect[3].y<<"  V_z  "<<vect[3].z<<endl;
	
	//cout<<"| "<<vect[2].x<<" "<<vect[2].y<<" "<<vect[2].z<<" 0 |"<<endl;
	//cout<<"| "<<vect[3].x<<" "<<vect[3].y<<" "<<vect[3].z<<" 0 |"<<endl;
	//cout<<"| "<<vect[1].x<<" "<<vect[1].y<<" "<<vect[1].z<<" 0 |"<<endl;
	//cout<<"| 0"<<" 0"<<" 0"<<" 1 |"<<endl;
	
	
	
	int index=0;
	int index2=0;
	int index_p=0;
	float u_img,v_img,m;
	float smin;	
	float tmin;
	int inter_obj=-1;
	float sphere_r,sphere_g,sphere_b,poly_r,poly_g,poly_b;
	float compare_s,compare_p;

	for(int x=0;x<scene[e].w;x++)
		for(int y=0;y<scene[e].h;y++)
		{
			
			float actual_x,actual_y,actual_z;
			int w=scene[e].w-1;
			int h=scene[e].h-1;
			
			//Compute image plane co-ordinated bounded in (-1,-1) to (1,1)
			u_img=float(x)*2/w-1;
			v_img=float(y)*2/h-1;
				
			float x1,y1,z1;
			x1=cam[f].p_x+(u_img*vect[2].x)+(v_img*vect[3].x)-(d*vect[1].x);
			y1=cam[f].p_y+(u_img*vect[2].y)+(v_img*vect[3].y)-(d*vect[1].y);
			z1=cam[f].p_z+(u_img*vect[2].z)+(v_img*vect[3].z)-(d*vect[1].z);
			
			m=sqrt(((x1-cam[f].p_x)*(x1-cam[f].p_x))+((y1-cam[f].p_y)*(y1-cam[f].p_y))+((z1-cam[f].p_z)*(z1-cam[f].p_z)));
		
			smin=cam[f].far;
			tmin=cam[f].far;
			
			index=1000;
			index_p=1000;
		
			//initialize to default
			c_pixel[x][y][0]=255;
			c_pixel[x][y][1]=255;
			c_pixel[x][y][2]=255;
			
			depth_pixel[x][y][0]=depth_pixel[x][y][1]=depth_pixel[x][y][2]=0;
			
			//Check for every object
			if(s_no!=0)//check if sphere objects exists in scene
			for(int i=0;i<s_no;i++)
			{
								
				//Create c vectore-Spehere centre to Ray start
				//Point- Point is vector
				float c_x,c_y,c_z;
				c_x=sph[i].X-cam[f].p_x;
				c_y=sph[i].Y-cam[f].p_y;
				c_z=sph[i].Z-cam[f].p_z;
				
				//Ray vector
				float v_x,v_y,v_z;
				
				v_x=(((x1-cam[f].p_x))/m)*c_x;
				v_y=(((y1-cam[f].p_y))/m)*c_y;
				v_z=(((z1-cam[f].p_z))/m)*c_z;
				sph[i].v=v_x+v_y+v_z;
				//Calculate v and  c
				//sph[i].v=sqrt(((v_x)*(v_x))+((v_y)*(v_y))+((v_z)*(v_z)));
				//calculate c
				sph[i].c=sqrt(((sph[i].X-cam[f].p_x)*(sph[i].X-cam[f].p_x))+((sph[i].Y-cam[f].p_y)*(sph[i].Y-cam[f].p_y))+((sph[i].Z-cam[f].p_z)*(sph[i].Z-cam[f].p_z)));
			
				//Calculate b			
				sph[i].b=sqrt((sph[i].c*sph[i].c)-(sph[i].v*sph[i].v));
				
			
				//Calculate d sqaure
				
				sph[i].dsquare=(sph[i].radius*sph[i].radius)-(sph[i].b*sph[i].b);
				
				
				if(sph[i].dsquare>=0)
				{
				
					
				
					//calculate s=(v-d)
					sph[i].s=sph[i].v-(sqrt(sph[i].dsquare));			
					
				
					if(sph[i].s<smin && sph[i].s>m)//check
					{
						
						smin=sph[i].s;								
						index=i;
						//cout<<"index "<<index<<endl;
										
					}
								
				}
				if(index==i)
				{
					inter_obj=index;

					/////1
					double I_r,I_g,I_b;
					double Ia=20;
					int colo=0;
					while((sph[index].usemtl).compare(material[colo].usemtl))
						colo++;
						
					//cout<<"Material "<<material[colo].usemtl<<" "<<material[colo].ks_r<<" "<<material[colo].ks_g<<" "<<material[colo].ks_b<<" "<<material[colo].kd_r<<" "<<material[colo].kd_g<<" "<<material[colo].kd_b<<" "<<material[colo].ka_r<<" "<<material[colo].ka_g<<" "<<material[colo].ka_b<<endl;


					//Ambient reflection
					I_r=Ia*material[colo].ka_r;
					I_g=Ia*material[colo].ka_g;
					I_b=Ia*material[colo].ka_b;
					
					//Diffuse Reflection
					
					//Ip*Kd*(N.L)
					
					
					float N_x,N_y,N_z,L_x,L_y,L_z;
					float N_mag;
					float L_mag;
					float diffuse,total_diffuse;
					float zero=0;
					double range=255;
					float V_x,V_y,V_z;
					float m2,m3;					
					float R_x,R_y,R_z;
					float V_dot_R;
					float specular;
					
					int shadow=0;
					float diff_r=0,diff_b=0,diff_g=0,spec_r=0,spec_g=0,spec_b=0;
					//find co-ordinates of	the point (x,y,z) on the sphere
					//PRP to (x1,y1,z1) plus s gives (x,y,z)

					actual_x=cam[f].p_x+(sph[index].s*(x1-cam[f].p_x))/m;
					actual_y=cam[f].p_y+(sph[index].s*(y1-cam[f].p_y))/m;
					actual_z=cam[f].p_z+(sph[index].s*(z1-cam[f].p_z))/m;


					//For Sphere N=(S-Pc)/|S-Pc|
					//L=P-S/|P-S|
					N_mag=sqrt(((actual_x-sph[index].X)*(actual_x-sph[index].X))+((actual_y-sph[index].Y)*(actual_y-sph[index].Y))+((actual_z-sph[index].Z)*(actual_z-sph[index].Z)));
					N_x=(actual_x-sph[index].X)/N_mag;
					N_y=(actual_y-sph[index].Y)/N_mag;
					N_z=(actual_z-sph[index].Z)/N_mag;//Surface Normal
					
					//calculate V=C-Q/|C-Q|
					
					m2=sqrt(((actual_x-cam[f].p_x)*(actual_x-cam[f].p_x))+((actual_y-cam[f].p_y)*(actual_y-cam[f].p_y))+((actual_z-cam[f].p_z)*(actual_z-cam[f].p_z)));
					V_x=(cam[f].p_x-actual_x)/m2;
					V_y=(cam[f].p_y-actual_y)/m2;
					V_z=(cam[f].p_z-actual_z)/m2;
					//cout<<" V_x"<<V_x<<" V_y "<<V_y<<" V_z "<< V_z<<endl;
					
					//For each light
					for(int ll=0;ll<l_no;ll++)
					{
						
					L_mag=sqrt(((lights[ll].x-actual_x)*(lights[ll].x-actual_x))+((lights[ll].y-actual_y)*(lights[ll].y-actual_y))+((lights[ll].z-actual_z)*(lights[ll].z-actual_z)));
					L_x=(lights[ll].x-actual_x)/L_mag;
					L_y=(lights[ll].y-actual_y)/L_mag;
					L_z=(lights[ll].z-actual_z)/L_mag;
			
					//Dot product of N and L
					//check dot product value and apply conditions
					/*
					 * //Dot product of N and L
					//check dot product value and apply conditions
					N_L=L_x*N_x+L_y*N_y+L_z*N_z;	
					if(N_L>=0){
					 */ 
					 
					 shadow=findShadow(actual_x,actual_y,actual_z,L_x,L_y,L_z,L_mag,inter_obj);
					//shadow=0;
					//if background color
					//intersect light with each sphere
					//PixelP0  Light(Lx,Ly,Lz)
				
					
					
				
					if(shadow)
					{
						//dnt add light
					}
					else
					{
					diffuse=(N_x*L_x)+(N_y*L_y)+(N_z*L_z);
					
					//if(diffuse<0)
						//diffuse=diffuse*(-1);
					
					diffuse=max(diffuse,zero);
					//cout<<"diffue "<<diffuse<<endl;
					
					I_r=I_r+(lights[ll].Ipx*material[colo].kd_r*diffuse);
					I_g=I_g+(lights[ll].Ipy*material[colo].kd_g*diffuse);
					I_b=I_b+(lights[ll].Ipz*material[colo].kd_b*diffuse);
					
					//Calculate Specular highlight
					
					//Calculate R=2(L.N)N-L
					R_x=(2*diffuse*N_x)-L_x;
					R_y=(2*diffuse*N_y)-L_y;
					R_z=(2*diffuse*N_z)-L_z;
					
					
					//Dot product of V and R
					V_dot_R=(V_x*R_x)+(V_y*R_y)+(V_z*R_z);
					V_dot_R=max(V_dot_R,zero);
					
					specular=pow(V_dot_R,material[colo].alpa);
					
					I_r=I_r+(lights[ll].Ipx*material[colo].ks_r*specular);
					I_g=I_g+(lights[ll].Ipy*material[colo].ks_g*specular);
					I_b=I_b+(lights[ll].Ipz*material[colo].ks_b*specular);
					
				
					}
				}

					I_r=min(I_r,range);					
					I_g=min(I_g,range);
					I_b=min(I_b,range);
					
					//set color for the surface of that sphere to the pixel***
					//c_pixel[x][y][0]=(int)I_r;
					//c_pixel[x][y][1]=(int)I_g;
					//c_pixel[x][y][2]=(int)I_b;
					
					sphere_r=I_r;
					sphere_g=I_g;
					sphere_b=I_b;
				}
			
		
			}
		
		
			float beta,gamma,tt;
			int index4=cam[f].far;
			
			
			if(f_no!=0)//Check for each polygonal face
			for(int jj=0;jj<f_no;jj++)
			{	
				float actual_x,actual_y,actual_z;
				int v0,v1,v2;
				v0=face[jj].a-1;
				v1=face[jj].b-1;
				v2=face[jj].c-1;
				//Get the three point of triangle
				vert[0].x=vertices[v0].x;
				vert[0].y=vertices[v0].y;
				vert[0].z=vertices[v0].z;
 
				vert[1].x=vertices[v1].x;
				vert[1].y=vertices[v1].y;
				vert[1].z=vertices[v1].z;
 
				vert[2].x=vertices[v2].x;
				vert[2].y=vertices[v2].y;
				vert[2].z=vertices[v2].z;
				
				
	
				//3 euations in 3 unknowns
				//L=>PRP
				//U=>x1-PRP/|x1-PRP|
				//
 
				float aa,bb,cc,dd,ee,ff,gg,hh,ii,pp,qq,rr;
				 //(B-A)===>A(a,b,c)
				aa=vert[3].x=vert[1].x-vert[0].x;
				dd=vert[3].y=vert[1].y-vert[0].y;
				gg=vert[3].z=vert[1].z-vert[0].z;
				//(C-A)==>B(d,e,f)
				bb=vert[4].x=vert[2].x-vert[0].x;
				ee=vert[4].y=vert[2].y-vert[0].y;
				hh=vert[4].z=vert[2].z-vert[0].z;
				//C==>-U(g,h,i)
				float Ux,Uy,Uz;
				float mag;
				mag=sqrt(((x1-cam[f].p_x)*(x1-cam[f].p_x))+((y1-cam[f].p_y)*(y1-cam[f].p_y))+((z1-cam[f].p_z)*(z1-cam[f].p_z)));
				cc=Ux=-(x1-cam[f].p_x)/mag;
				ff=Uy=-(y1-cam[f].p_y)/mag;
				ii=Uz=-(z1-cam[f].p_z)/mag;

				//D==>L-A(p,q,r)
				float Dx,Dy,Dz;
				pp=Dx=cam[f].p_x-vert[0].x;
				qq=Dy=cam[f].p_y-vert[0].y;
				rr=Dz=cam[f].p_z-vert[0].z;

				//find beta, gamma and t

				float dtr=aa*ee*ii+bb*ff*gg+cc*dd*hh-(cc*ee*gg+aa*ff*hh+dd*bb*ii);	
				float aaa=ee*ii-ff*hh;
				float bbb=cc*hh-bb*ii;
				float ccc=bb*ff-cc*ee;
				float ddd=ff*gg-dd*ii;
				float eee=aa*ii-cc*gg;
				float fff=cc*dd-aa*ff;
				float ggg=dd*hh-ee*gg;
				float hhh=bb*gg-aa*hh;
				float iii=aa*ee-bb*dd;
				if(dtr!=0){
					beta=(pp*aaa+qq*bbb+rr*ccc)/dtr;
					gamma=(pp*ddd+qq*eee+rr*fff)/dtr;
					tt=(pp*ggg+qq*hhh+rr*iii)/dtr;
				}
				else{
					beta=0;
					gamma=0;
					tt=0;
		
				}

				if(beta>=0 && gamma>=0 && (beta+gamma)<=1)
				{
					if(tt<tmin){
						tmin=tt;
						index_p=jj;
					}
					
			
				}
				if(index_p==jj)
				{
					inter_obj=index_p;
					
					vect[5].x=vect[3].x;
					vect[5].y=vect[3].y;
					vect[5].z=vect[3].z;
					
					vect[6].x=vect[4].x;
					vect[6].y=vect[4].y;
					vect[6].z=vect[4].z;
					
					double I_r,I_g,I_b;
					double Ia=20;
					int  count;
					int colo=0;
					if(!face[1].usemtl.compare("NA")){
						count=0;
					}
					else
						count=index_p;
						
					while((face[count].usemtl).compare(material[colo].usemtl))//***change it for each face
						colo++;
					//cout<<face[0].usemtl<<endl;
					//cout<<"Material "<<material[colo].usemtl<<" "<<material[colo].ks_r<<" "<<material[colo].ks_g<<" "<<material[colo].ks_b<<" "<<material[colo].kd_r<<" "<<material[colo].kd_g<<" "<<material[colo].kd_b<<" "<<material[colo].ka_r<<" "<<material[colo].ka_g<<" "<<material[colo].ka_b<<endl;
					//Ambient reflection
					I_r=Ia*material[colo].ka_r;
					I_g=Ia*material[colo].ka_g;
					I_b=Ia*material[colo].ka_b;
					
					//Diffuse Reflection
					
					//Ip*Kd*(N.L)
					
					
					float N_x,N_y,N_z,L_x,L_y,L_z;
					float N_mag;
					float L_mag;
					float diffuse,total_diffuse;
					float zero=0;
					double range=255;
					float V_x,V_y,V_z;
					float m2,m3;					
					float R_x,R_y,R_z;
					float V_dot_R;
					float specular;
					float diff_r=0,diff_b=0,diff_g=0,spec_r=0,spec_g=0,spec_b=0;
					int shadow=0;
					//find co-ordinates of	the point (x,y,z) on the sphere
					//PRP to (x1,y1,z1) plus s gives (x,y,z)

					actual_x=cam[f].p_x+(tmin*(x1-cam[f].p_x))/m;
					actual_y=cam[f].p_y+(tmin*(y1-cam[f].p_y))/m;
					actual_z=cam[f].p_z+(tmin*(z1-cam[f].p_z))/m;

					//For polygon N=E1xE2/|E1xE2|
					
					
					/*
					N_mag=sqrt(((actual_x-sph[index].X)*(actual_x-sph[index].X))+((actual_y-sph[index].Y)*(actual_y-sph[index].Y))+((actual_z-sph[index].Z)*(actual_z-sph[index].Z)));
					N_x=(actual_x-sph[index].X)/N_mag;
					N_y=(actual_y-sph[index].Y)/N_mag;
					N_z=(actual_z-sph[index].Z)/N_mag;//Surface Normal
					* */
					
					
					////
					
					
					//Vert[3]=E1 Vert[4]=E2
					N_x = vert[3].y*vert[4].z - vert[4].y*vert[3].z;	
					N_y =  vert[4].x*vert[3].z-vert[3].x*vert[4].z;
					N_z = vert[3].x*vert[4].y - vert[3].y*vert[4].x; 
					N_mag=sqrt((N_x*N_x)+(N_y*N_y)+(N_z*N_z));
					N_x=N_x/N_mag;
					N_y=N_y/N_mag;
					N_z=N_z/N_mag;
					
					//calculate V=C-Q/|C-Q|
					
					m2=sqrt(((actual_x-cam[f].p_x)*(actual_x-cam[f].p_x))+((actual_y-cam[f].p_y)*(actual_y-cam[f].p_y))+((actual_z-cam[f].p_z)*(actual_z-cam[f].p_z)));
					V_x=(cam[f].p_x-actual_x)/m2;
					V_y=(cam[f].p_y-actual_y)/m2;
					V_z=(cam[f].p_z-actual_z)/m2;
					//cout<<" V_x"<<V_x<<" V_y "<<V_y<<" V_z "<< V_z<<endl;
					
					//For each light
					for(int ll=0;ll<l_no;ll++)
					{
					float N_L;		
					L_mag=sqrt(((lights[ll].x-actual_x)*(lights[ll].x-actual_x))+((lights[ll].y-actual_y)*(lights[ll].y-actual_y))+((lights[ll].z-actual_z)*(lights[ll].z-actual_z)));
					L_x=(lights[ll].x-actual_x)/L_mag;
					L_y=(lights[ll].y-actual_y)/L_mag;
					L_z=(lights[ll].z-actual_z)/L_mag;
			
					//Dot product of N and L
					//check dot product value and apply conditions
					//Dot product of N and L
					//check dot product value and apply conditions
					N_L=L_x*N_x+L_y*N_y+L_z*N_z;	
					//if(N_L>=0){
					//if the pixel is in shadow
					
					
					shadow=findShadow(actual_x,actual_y,actual_z,L_x,L_y,L_z,L_mag,inter_obj);
					//shadow=0;
					//if background color
					//intersect light with each sphere
					//PixelP0  Light(Lx,Ly,Lz)
				
					
					
				
					if(shadow)
					{
						//dnt add light
					}
					else
					{
					diffuse=(N_x*L_x)+(N_y*L_y)+(N_z*L_z);
					
					//if(diffuse<0)
						//diffuse=diffuse*(-1);
					
					diffuse=max(diffuse,zero);
					//cout<<"diffue "<<diffuse<<endl;
					
					I_r=I_r+(lights[ll].Ipx*material[colo].kd_r*diffuse);
					I_g=I_g+(lights[ll].Ipy*material[colo].kd_g*diffuse);
					I_b=I_b+(lights[ll].Ipz*material[colo].kd_b*diffuse);
					
					//Calculate Specular highlight
					
					//Calculate R=2(L.N)N-L
					R_x=(2*diffuse*N_x)-L_x;
					R_y=(2*diffuse*N_y)-L_y;
					R_z=(2*diffuse*N_z)-L_z;
					
					
					//Dot product of V and R
					V_dot_R=(V_x*R_x)+(V_y*R_y)+(V_z*R_z);
					V_dot_R=max(V_dot_R,zero);
					
					specular=pow(V_dot_R,material[colo].alpa);
					
					I_r=I_r+(lights[ll].Ipx*material[colo].ks_r*specular);
					I_g=I_g+(lights[ll].Ipy*material[colo].ks_g*specular);
					I_b=I_b+(lights[ll].Ipz*material[colo].ks_b*specular);
					
				
					}


					I_r=min(I_r,range);					
					I_g=min(I_g,range);
					I_b=min(I_b,range);
					
					poly_r=I_r;
					poly_g=I_g;
					poly_b=I_b;
				}
					
				///2
				}
			
			}

			//Intersection with nearest object done
			string object;
			float near_t;
			if(inter_obj!=-1)
			if(smin>tmin){
			c_pixel[x][y][0]=(int)poly_r;
			c_pixel[x][y][1]=(int)poly_g;
			c_pixel[x][y][2]=(int)poly_b;
			//calulate depth_value
					float p=255*(tmin-m);//index4 changed to tmin
					float dist=p/(cam[f].far-m);
					int depth=255-(min((float)255,dist));
					depth_pixel[x][y][0]=depth_pixel[x][y][1]=depth_pixel[x][y][2]=depth;
					near_t=tmin;
					object="poly";
			}
			else if (smin<tmin){
			c_pixel[x][y][0]=(int)sphere_r;
			c_pixel[x][y][1]=(int)sphere_g;
			c_pixel[x][y][2]=(int)sphere_b;
			
			//calulate depth_value
			float p=255*(smin-m);
			float dist=p/(cam[f].far-m);
			int depth=255-(min((float)255,dist));
			depth_pixel[x][y][0]=depth_pixel[x][y][1]=depth_pixel[x][y][2]=depth;
			near_t=smin;
			object="sphere";
			////3
			}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//near_t,object,inter_obj
			/*
			float actual_x,actual_y,actual_z;
			actual_x=cam[f].p_x+(near_t*(x1-cam[f].p_x))/m;
			actual_y=cam[f].p_y+(near_t*(y1-cam[f].p_y))/m;
			actual_z=cam[f].p_z+(near_t*(z1-cam[f].p_z))/m;
					
				
			
			double I_r=0,I_g=0,I_b=0;
			double Ia=20;
			int colo=0;
			if(!object.compare("sphere"))
				while((sph[inter_obj].usemtl).compare(material[colo].usemtl))
							colo++;
			if(!object.compare("poly"))
				while((face[0].usemtl).compare(material[colo].usemtl))//**change 0 to actaul
							colo++;	
			
			if(!object.compare("poly"))
				cout<<"Material "<<material[colo].usemtl<<" "<<material[colo].ks_r<<" "<<material[colo].ks_g<<" "<<material[colo].ks_b<<" "<<material[colo].kd_r<<" "<<material[colo].kd_g<<" "<<material[colo].kd_b<<" "<<material[colo].ka_r<<" "<<material[colo].ka_g<<" "<<material[colo].ka_b<<endl;
			//cout<<"Material "<<colo<<endl;;
	
			//cout<<"Inetersection"<<object<<endl;
			if(colo!=-1){
			//Ambient reflection
			
			I_r=Ia*material[colo].ka_r;
			I_g=Ia*material[colo].ka_g;
			I_b=Ia*material[colo].ka_b;
			}		
			//Diffuse Reflection
					
			//Ip*Kd*(N.L)
					
					
			float N_x,N_y,N_z,L_x,L_y,L_z;
			float N_mag;
			float L_mag;
			float diffuse,total_diffuse;
			float zero=0;
			double range=255;
			float V_x,V_y,V_z;
			float m2,m3;					
			float R_x,R_y,R_z;
			float L_N;
			float V_dot_R;
			float specular;

			int shadow=0;
			
			float diff_r=0,diff_b=0,diff_g=0,spec_r=0,spec_g=0,spec_b=0;
			if(inter_obj!=-1)
			if(!object.compare("sphere")){
					


					N_mag=sqrt(((actual_x-sph[inter_obj].X)*(actual_x-sph[inter_obj].X))+((actual_y-sph[inter_obj].Y)*(actual_y-sph[inter_obj].Y))+((actual_z-sph[inter_obj].Z)*(actual_z-sph[inter_obj].Z)));
					N_x=(actual_x-sph[inter_obj].X)/N_mag;
					N_y=(actual_y-sph[inter_obj].Y)/N_mag;
					N_z=(actual_z-sph[inter_obj].Z)/N_mag;//Surface Normal
					
					//calculate V=C-Q/|C-Q|
					
					m2=sqrt(((actual_x-cam[f].p_x)*(actual_x-cam[f].p_x))+((actual_y-cam[f].p_y)*(actual_y-cam[f].p_y))+((actual_z-cam[f].p_z)*(actual_z-cam[f].p_z)));
					V_x=(cam[f].p_x-actual_x)/m2;
					V_y=(cam[f].p_y-actual_y)/m2;
					V_z=(cam[f].p_z-actual_z)/m2;
					//cout<<" V_x"<<V_x<<" V_y "<<V_y<<" V_z "<< V_z<<endl;
					
					bool isInShadow=false;
					//index2=1000;
					//float smin2=cam2[f].far;
					
					//For each light
					for(int ll=0;ll<l_no;ll++)
					{
					float N_L;	
					L_mag=sqrt(((lights[ll].x-actual_x)*(lights[ll].x-actual_x))+((lights[ll].y-actual_y)*(lights[ll].y-actual_y))+((lights[ll].z-actual_z)*(lights[ll].z-actual_z)));
					L_x=(lights[ll].x-actual_x)/L_mag;
					L_y=(lights[ll].y-actual_y)/L_mag;
					L_z=(lights[ll].z-actual_z)/L_mag;
					
					
					//Dot product of N and L
					//check dot product value and apply conditions
					N_L=L_x*N_x+L_y*N_y+L_z*N_z;	
					if(N_L>=0){
					//if the pixel is in shadow
					shadow=findShadow(actual_x,actual_y,actual_z,L_x,L_y,L_z,L_mag,inter_obj);
					//shadow=0;
					//if background color
					//intersect light with each sphere
					//PixelP0  Light(Lx,Ly,Lz)
				
					
					
				
					if(shadow)
					{
						//dnt add light
					}
					else
					{
					diffuse=(N_x*L_x)+(N_y*L_y)+(N_z*L_z);
					
					//if(diffuse<0)
						//diffuse=diffuse*(-1);
					
					diffuse=max(diffuse,zero);
					//cout<<"diffue "<<diffuse<<endl;
					
					I_r=I_r+(lights[ll].Ipx*material[colo].kd_r*diffuse);
					I_g=I_g+(lights[ll].Ipy*material[colo].kd_g*diffuse);
					I_b=I_b+(lights[ll].Ipz*material[colo].kd_b*diffuse);
					
					//Calculate Specular highlight
					
					//Calculate R=2(L.N)N-L
					R_x=(2*diffuse*N_x)-L_x;
					R_y=(2*diffuse*N_y)-L_y;
					R_z=(2*diffuse*N_z)-L_z;
					
					
					//Dot product of V and R
					V_dot_R=(V_x*R_x)+(V_y*R_y)+(V_z*R_z);
					V_dot_R=max(V_dot_R,zero);
					
					specular=pow(V_dot_R,material[colo].alpa);
					
					I_r=I_r+(lights[ll].Ipx*material[colo].ks_r*specular);
					I_g=I_g+(lights[ll].Ipy*material[colo].ks_g*specular);
					I_b=I_b+(lights[ll].Ipz*material[colo].ks_b*specular);
					}
				}
					}


					I_r=min(I_r,range);					
					I_g=min(I_g,range);
					I_b=min(I_b,range);
					
					//set color for the surface of that sphere to the pixel
					c_pixel[x][y][0]=(int)I_r;
					c_pixel[x][y][1]=(int)I_g;
					c_pixel[x][y][2]=(int)I_b;
					
					
			}
			else if(!object.compare("poly")){
					
				
				/*
					//N=E1XE2/|E1XE2|
					
					//E1=Vect[5]  E2=Vect[6]
					//U-->Vect[2]
					N_x = vect[5].y*vect[6].z - vect[6].y*vect[5].z;	
					N_y=  vect[6].x*vect[5].z-vect[5].x*vect[6].z;
					N_z  = vect[5].x*vect[6].y - vect[5].y*vect[6].x; 
									
					N_mag=sqrt((N_x*N_x)+(N_y*N_y)+(N_z*N_z));
					if(N_mag!=0){
					N_x=(N_x/N_mag);
					N_y=(N_y/N_mag);
					N_z=(N_z/N_mag);
					}
				*/
				
				
				
				/*
				//Vert[3]=E1 Vert[4]=E2
					N_x = vert[3].y*vert[4].z - vert[4].y*vert[3].z;	
					N_y =  vert[4].x*vert[3].z-vert[3].x*vert[4].z;
					N_z = vert[3].x*vert[4].y - vert[3].y*vert[4].x; 
					N_mag=sqrt((N_x*N_x)+(N_y*N_y)+(N_z*N_z));
					N_x=N_x/N_mag;
					N_y=N_y/N_mag;
					N_z=N_z/N_mag;
					//calculate V=C-Q/|C-Q|
					
					m2=sqrt(((actual_x-cam[f].p_x)*(actual_x-cam[f].p_x))+((actual_y-cam[f].p_y)*(actual_y-cam[f].p_y))+((actual_z-cam[f].p_z)*(actual_z-cam[f].p_z)));
					V_x=(cam[f].p_x-actual_x)/m2;
					V_y=(cam[f].p_y-actual_y)/m2;
					V_z=(cam[f].p_z-actual_z)/m2;
					//cout<<" V_x"<<V_x<<" V_y "<<V_y<<" V_z "<< V_z<<endl;
					
					bool isInShadow=false;
					//index2=1000;
					//float smin2=cam2[f].far;
					
					//For each light
					for(int ll=0;ll<l_no;ll++)
					{
					float N_L;	
					
					//Ray
					L_mag=sqrt(((lights[ll].x-actual_x)*(lights[ll].x-actual_x))+((lights[ll].y-actual_y)*(lights[ll].y-actual_y))+((lights[ll].z-actual_z)*(lights[ll].z-actual_z)));
					L_x=(lights[ll].x-actual_x)/L_mag;
					L_y=(lights[ll].y-actual_y)/L_mag;
					L_z=(lights[ll].z-actual_z)/L_mag;
					
					
					//Dot product of N and L
					//check dot product value and apply conditions
					N_L=L_x*N_x+L_y*N_y+L_z*N_z;	
					if(N_L>=0){
					//if the pixel is in shadow
					
					
					//shadow=findShadow(actual_x,actual_y,actual_z,L_x,L_y,L_z,L_mag,inter_obj);
					shadow=0;
					//if background color
					//intersect light with each sphere
					//PixelP0  Light(Lx,Ly,Lz)
				
					
					
				
					if(shadow)
					{
						//dnt add light
					}
					else
					{
					diffuse=(N_x*L_x)+(N_y*L_y)+(N_z*L_z);
					
					//if(diffuse<0)
						//diffuse=diffuse*(-1);
					
					diffuse=max(diffuse,zero);
					//cout<<"diffue "<<diffuse<<endl;
					
					I_r=I_r+(lights[ll].Ipx*material[colo].kd_r*diffuse);
					I_g=I_g+(lights[ll].Ipy*material[colo].kd_g*diffuse);
					I_b=I_b+(lights[ll].Ipz*material[colo].kd_b*diffuse);
					
					//Calculate Specular highlight
					
					//Calculate R=2(L.N)N-L
					R_x=(2*diffuse*N_x)-L_x;
					R_y=(2*diffuse*N_y)-L_y;
					R_z=(2*diffuse*N_z)-L_z;
					
					
					//Dot product of V and R
					V_dot_R=(V_x*R_x)+(V_y*R_y)+(V_z*R_z);
					V_dot_R=max(V_dot_R,zero);
					
					specular=pow(V_dot_R,material[colo].alpa);
					
					I_r=I_r+(lights[ll].Ipx*material[colo].ks_r*specular);
					I_g=I_g+(lights[ll].Ipy*material[colo].ks_g*specular);
					I_b=I_b+(lights[ll].Ipz*material[colo].ks_b*specular);
					}
				}
					}


					I_r=min(I_r,range);					
					I_g=min(I_g,range);
					I_b=min(I_b,range);
					
					//set color for the surface of that sphere to the pixel
					c_pixel[x][y][0]=(int)I_r;
					c_pixel[x][y][1]=(int)I_g;
					c_pixel[x][y][2]=(int)I_b;
					
			
			}
			*/

			///////////////////////////////////
		}	
	
		//Output the color_values and depth_values to save in ppm
		img_w=scene[e].w;
		img_h=scene[e].h;
		saveToPPM(e,f);
}

int main(int argc,char* argv[])
{
	fileParsing(argv);
	
	for(int a=0;a<r_no;a++)
	{
		
		//for each scene
		for(int b=0;b<c_no;b++)
		{
			//for each camera
			trace(a,b);
			
		}
	}
	return 0;
}
