#include <iostream>
#include <sstream>
#include<vector>
#include <fstream>
#include <cstring>
#include <string>
#include<math.h>
#include<algorithm>

using namespace std;


using std::cout;
using std::endl;
using std::istringstream;
using std::ifstream;

//Data structure for Ray,Sphere,Camera
struct Camera{
string camera_name;
float p_x,p_y,p_z,v_x,v_y,v_z;
float near,far;
};

struct Ray{
string scene_cast;
int w,h,recurr;
};

struct Sphere{
string sphere_name;
float X,Y,Z;
int R,G,B;
float radius;
float v,c,b,dsquare,s;
};


Camera *cam=new Camera[1000];
Ray *scene=new Ray[1000];
Sphere *sph=new Sphere[1000];
int c_no,s_no,r_no;
int img_w=0,img_h=0;
int c_pixel[5000][5000][3];
int depth_pixel[5000][5000][3];



void fileParsing(char *argv[])
{
	c_no=s_no=r_no=0;//No of Camera,Sphere and Ray
	int next;
	int int_obj;
	float float_obj;
	

	std::vector<std::string> parsed_token;
	std::vector<std::string> parsed_token1;
	std::string str;
	
	
	// create a file-reading object
	ifstream fin;
	fin.open(argv[1],std::ifstream::in);//open commands01.txt
	
	if ( !fin.is_open() )
      cout<<"Could not open file\n";
	else
	{
		
		while (fin >> str)
			parsed_token.push_back(str);
		for (int i = 0; i < parsed_token.size(); ++i)
		{
			
			if(parsed_token[i]=="s")
			{
				istringstream s3;
				s3.str(parsed_token[i+1]);
				s3>>sph[s_no].sphere_name;
				
				
				if ( ! (istringstream(parsed_token[i+2]) >> float_obj) ) float_obj = 0;
					sph[s_no].X=float_obj;
				if ( ! (istringstream(parsed_token[i+3]) >> float_obj) ) float_obj = 0;
					sph[s_no].Y=float_obj;
				if ( ! (istringstream(parsed_token[i+4]) >> float_obj) ) float_obj = 0;
					sph[s_no].Z=float_obj;
			
				if ( ! (istringstream(parsed_token[i+5]) >> float_obj) ) float_obj = 0;
					sph[s_no].radius=float_obj;

				if ( ! (istringstream(parsed_token[i+6]) >> int_obj) ) int_obj = 0;
					sph[s_no].R=int_obj;
				if ( ! (istringstream(parsed_token[i+7]) >> int_obj) ) int_obj = 0;
					sph[s_no].G=int_obj;
				if ( ! (istringstream(parsed_token[i+8]) >> int_obj) ) int_obj = 0;
					sph[s_no].B=int_obj;
			
				next=i+8;
				s_no++;
				i=next;
			
			}
		}
		
  }
	fin.close();//commands01.txt file reading finished
	
	//creating file reading object
	ifstream obj;
	obj.open(argv[2],std::ifstream::in); //filename object file
                                   
    if ( !obj.is_open() )
      cout<<"Could not open file\n";
	else
	{
	
		while (obj >> str)
			parsed_token1.push_back(str);
		for (int i = 0; i < parsed_token1.size(); ++i)
		{
			if(parsed_token1[i]=="c")
			{
				
				istringstream s1;
				s1.str(parsed_token1[i+1]);
				s1>>cam[c_no].camera_name;
				
				
				if ( ! (istringstream(parsed_token1[i+2]) >> float_obj) ) float_obj = 0;
					cam[c_no].p_x=float_obj;
				if ( ! (istringstream(parsed_token1[i+3]) >> float_obj) ) float_obj = 0;
					cam[c_no].p_y=float_obj;
				if ( ! (istringstream(parsed_token1[i+4]) >> float_obj) ) float_obj = 0;
					cam[c_no].p_z=float_obj;
			
				if ( ! (istringstream(parsed_token1[i+5]) >> float_obj) ) float_obj = 0;
					cam[c_no].v_x=float_obj;
				if ( ! (istringstream(parsed_token1[i+6]) >> float_obj) ) float_obj = 0;
					cam[c_no].v_y=float_obj;
				if ( ! (istringstream(parsed_token1[i+7]) >> float_obj) ) float_obj = 0;
					cam[c_no].v_z=float_obj;
				
				if ( ! (istringstream(parsed_token1[i+8]) >> float_obj) ) float_obj = 0;
					cam[c_no].near=float_obj;
				if ( ! (istringstream(parsed_token1[i+9]) >> float_obj) ) float_obj = 0;
					cam[c_no].far=float_obj;
			
				
				c_no++;
				next=i+9;
				i=next;
			
			}
			else if(parsed_token1[i]=="r")
			{
				istringstream s2;
				s2.str(parsed_token1[i+1]);
				s2>>scene[r_no].scene_cast;
				
				if ( ! (istringstream(parsed_token1[i+2]) >> int_obj) ) int_obj = 0;
					scene[r_no].w=int_obj;
				if ( ! (istringstream(parsed_token1[i+3]) >> int_obj) ) int_obj = 0;
					scene[r_no].h=int_obj;
				if ( ! (istringstream(parsed_token1[i+4]) >> int_obj) ) int_obj = 0;
					scene[r_no].recurr=int_obj;
			
				
				r_no++;
				next=i+4;
				i=next;
			
			}
		
		}
	}
	obj.close();//finished parsing obj file
	
	
}

void saveToPPM(int e,int f){

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
	
	//Depth image
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


void trace(int e,int f)
{
	
	int index=0;//track which sphere is hit at the pixel by the ray***
	float u_img,v_img,m;
	float smin;



	for(int x=0;x<scene[e].w;x++)//made <= to get last bound(1,1)***
		for(int y=0;y<scene[e].h;y++)
		{
			
			
			int w=scene[e].w-1;//check***
			int h=scene[e].h-1;//check***
			
			//Compute image plane co-ordinated bounded in (-1,-1) to (1,1)
			u_img=float(x)*2/w-1;//or made 2/(img_width-1) to get (1,1) bound//check***
			v_img=float(y)*2/h-1;//check***
			float w_img=-cam[f].near;
			
			//Caculate m
			m=sqrt(u_img*u_img+v_img*v_img+w_img*w_img);
			
			smin=cam[f].far;
			
			
			index=1000;//***
			
			//initialize to default
			c_pixel[x][y][0]=0;
			c_pixel[x][y][1]=0;
			c_pixel[x][y][2]=0;
			
			depth_pixel[x][y][0]=depth_pixel[x][y][1]=depth_pixel[x][y][2]=0;
			
			for(int i=0;i<s_no;i++)
			{
						
				//computing V-> dot product of Centre and image plane
				
				sph[i].v=(sph[i].X*(u_img/m))+(sph[i].Y*(v_img/m))+(sph[i].Z*(w_img/m));
				
			
	
				//Calculate C=>distance between centre of origin and prp
				
				sph[i].c=sqrt(((sph[i].X-cam[f].p_x)*(sph[i].X-cam[f].p_x))+((sph[i].Y-cam[f].p_y)*(sph[i].Y-cam[f].p_y))+((sph[i].Z-cam[f].p_z)*(sph[i].Z-cam[f].p_z)));
				
				
			
				//Calculate b			
				sph[i].b=sqrt((sph[i].c*sph[i].c)-(sph[i].v*sph[i].v));
				
			
				//Calculate d sqaure
				
				sph[i].dsquare=(sph[i].radius*sph[i].radius)-(sph[i].b*sph[i].b);
				
				
				if(sph[i].dsquare>0)
				{
				
					
				
					//calculate s=(v-d)
					sph[i].s=sph[i].v-(sqrt(sph[i].dsquare));			
					
				
					if(sph[i].s<smin && sph[i].s>m)
					{
						
						smin=sph[i].s;								
						index=i;
										
					}
								
				}
				if(index==i)
				{
					//set color for the surface of that sphere to the pixel
					c_pixel[x][y][0]=sph[index].R;
					c_pixel[x][y][1]=sph[index].G;
					c_pixel[x][y][2]=sph[index].B;
					
					//calulate depth_value
					float p=255*(smin-m);
					float dist=p/(cam[f].far-m);
					int depth=255-(min((float)255,dist));
					depth_pixel[x][y][0]=depth_pixel[x][y][1]=depth_pixel[x][y][2]=depth;
				}
			
			}
		
		
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
