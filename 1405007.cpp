#include <iostream>
#include <bits/stdc++.h>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include "bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))
#define epsilon (1.0e-6)

class homogeneous_point
{
public:
    double x, y, z, w;

    // set the three coordinates, set w to 1
    homogeneous_point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1;
    }

    /*
    default constructor. does nothing. allows declarations like below:
        matrix m;
    therefore, usage is dangerous
    */
    homogeneous_point() {
    }

    // constructs a homogeneous point with given coordinates. forces w to be 1.0
    // if w is zero, raises error
    homogeneous_point(double x, double y, double z, double w)
    {
        assert (w != 0);
        this->x = x/w;
        this->y = y/w;
        this->z = z/w;
        this->w = 1;
    }

    // adds two points. returns a point forcing w to be 1.0
    homogeneous_point operator+ (const homogeneous_point& point)
    {
        double x = this->x + point.x;
        double y = this->y + point.y;
        double z = this->z + point.z;
        double w = this->w + point.w;
        homogeneous_point p(x, y, z, w);
        return p;
    }

    // subtracts one point from another. returns a point forcing w to be 1.0
    homogeneous_point operator- (const homogeneous_point& point)
    {
        double x = this->x - point.x;
        double y = this->y - point.y;
        double z = this->z - point.z;
        double w = this->w - point.w;
        homogeneous_point p(x, y, z, w);
        return p;
    }

    bool operator== (const homogeneous_point& point)
    {
        return (this->x==point.x && this->y==point.y && this->z==point.z);
    }

    // Print the coordinates of a point. exists for testing purpose.
    void print()
    {
        cout << "Point: " << endl;
        cout << x << " " << y << " " << z << " " << w << endl;
    }

};


class Vector
{
public:
    double x, y, z;

    // constructs a vector with given components
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // keeps the direction same. recalculates the vector to be unit.
    void normalize()
    {
        double r = sqrt(x*x + y*y + z*z);
        x = x / r;
        y = y / r;
        z = z / r;
    }

    // add two vectors
    Vector operator+(const Vector& v)
    {
        Vector v1(x+v.x, y+v.y, z+v.z);
        return v1;
    }

    // subtract one vector from another
    Vector operator-(const Vector& v)
    {
        Vector v1(x-v.x, y-v.y, z-v.z);
        return v1;
    }

    // scale a vector with a given coefficient
    Vector operator* (double m)
    {
        Vector v(x*m, y*m, z*m);
        return v;
    }

    // get the dot product of two vectors
    static double dot(Vector a, Vector b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    // get the cross product of two vectors
    static Vector cross(Vector a, Vector b)
    {
        Vector v(a.y*b.z - a.z*b.y, b.x*a.z - b.z*a.x, a.x*b.y - a.y*b.x);
        return v;
    }

    // Given vector a (axix of rotation), x (the vector to be rotated), angle of rotation
    static Vector R(Vector x, Vector a, double angle)
    {
        return (x*cos(angle*2*pi/360)+a*((1-cos(angle*2*pi/360))*dot(a,x)))+cross(a,x)*sin(angle*2*pi/360);
    }

    // print a vector. only for testing purposes.
    void print ()
    {
        cout << "Vector" << endl;
        cout << x << " " << y << " " << z << endl;
    }
};


/*
The matrices are forced to be 4x4. This is because in this assignment, we will deal with points in triangles.
Maximum # of points that we will deal with at once is 3. And all the standard matrices are 4x4 (i.e. scale, translation, rotation etc.)
*/
class matrix
{
public:
    double values[4][4];
    int num_rows, num_cols;

    // only set the number of rows and cols
    matrix(int rows, int cols)
    {
        assert (rows <= 4 && cols <= 4);
        num_rows = rows;
        num_cols = cols;
    }

    // prepare an nxn square matrix
    matrix(int n)
    {
        assert (n <= 4);
        num_rows = num_cols = n;
    }

    // prepare and return an identity matrix of size nxn
    static matrix make_identity(int n)
    {
        assert (n <= 4);
        matrix m(n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                    m.values[i][j] = 1;
                else
                    m.values[i][j] = 0;
            }
        }
        return m;
    }

    // print the matrix. exists for testing purposes
    void print()
    {
        cout << "Matrix:" << endl;
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                cout << values[i][j] << "\t";
            }
            cout << endl;
        }
    }

    // add the two matrices. Raise error if dimension mismatches
    matrix operator+ (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] + m.values[i][j];
            }
        }
        return m1;
    }

    // subtract a matrix from another. raise error if dimension mismatches
    matrix operator- (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] - m.values[i][j];
            }
        }
        return m1;
    }

    // multiply two matrices. allows statements like m1 = m2 * m3; raises error is dimension mismatches
    matrix operator* (const matrix& m)
    {
        assert (this->num_cols == m.num_rows);
        matrix m1(this->num_rows, m.num_cols);

        for (int i = 0; i < m1.num_rows; i++) {
            for (int j = 0; j < m1.num_cols; j++) {
                double val = 0;
                for (int k = 0; k < this->num_cols; k++) {
                    val += this->values[i][k] * m.values[k][j];
                }
                m1.values[i][j] = val;
            }
        }
        return m1;
    }

    // multiply a matrix with a constant
    matrix operator* (double m)
    {
        matrix m1(this->num_rows, this->num_cols);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m1.values[i][j] = m * this->values[i][j];
            }
        }
        return m1;
    }

    // multiply a 4x4 matrix with a homogeneous point and return the resulting point.
    // usage: homogeneous_point p = m * p1;
    // here, m is a 4x4 matrix, intended to be the transformation matrix
    // p1 is the point on which the transformation is being made
    // p is the resulting homogeneous point
    homogeneous_point operator* (const homogeneous_point& p)
    {
        assert (this->num_rows == this->num_cols && this->num_rows == 4);

        matrix m(4, 1);
        m.values[0][0] = p.x;
        m.values[1][0] = p.y;
        m.values[2][0] = p.z;
        m.values[3][0] = p.w;

        matrix m1 = (*this)*m;
        homogeneous_point p1(m1.values[0][0], m1.values[1][0], m1.values[2][0], m1.values[3][0]);
        return p1;
    }

    // return the transpose of a matrix
    matrix transpose()
    {
        matrix m(num_cols, num_rows);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m.values[j][i] = values[i][j];
            }
        }
        return m;
    }

    //take tx, ty, tz and return a translation matrix
    static matrix translate(double tx, double ty, double tz)
    {
        matrix m(4);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                //m.values[j][i] = values[i][j];
                if(i==j)
                {
                    m.values[i][j]=1;
                }
                else{
                    m.values[i][j]=0;
                }
            }
        }
        m.values[0][3]=tx;
        m.values[1][3]=ty;
        m.values[2][3]=tz;

        return m;
    }

    //Given sx, sy, sz . return scale matrix
    static matrix scale(double sx, double sy, double sz)
    {
        matrix m(4);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                //m.values[j][i] = values[i][j];
                if(i==j)
                {
                    m.values[i][j]=1;
                }
                else{
                    m.values[i][j]=0;
                }
            }
        }
        m.values[0][0]=sx;
        m.values[1][1]=sy;
        m.values[2][2]=sz;

        return m;

    }

    //Zero matrix
    static matrix Zero()
    {
        matrix m(4);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m.values[i][j]=0;
            }
        }
        return m;


    }

};

/*
A simple class to hold the color components, r, g, b of a certain shade.
*/
class color {
public:
    double r, g, b;
    color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    color() {
    }
};

homogeneous_point Line_Plane_Intersect(homogeneous_point a, homogeneous_point b, double z_val)
{
    double t=(z_val-a.z)/(b.z-a.z);

    double rx=a.x+t*(b.x-a.x);
    double ry=a.y+t*(b.y-a.y);
    double rz=a.z+t*(b.z-a.z);

    homogeneous_point r(rx,ry,rz);

    return r;
}

homogeneous_point Line_Plane_Intersect1(homogeneous_point a, homogeneous_point b, double y_val)
{

    double t=(y_val-a.y)/(b.y-a.y);

    double rx=a.x+t*(b.x-a.x);
    double ry=a.y+t*(b.y-a.y);
    double rz=a.z+t*(b.z-a.z);

    homogeneous_point r(rx,ry,rz);

    return r;
}

homogeneous_point max_Y(homogeneous_point a, homogeneous_point b, homogeneous_point c)
{
    if(a.y>=b.y && a.y>=c.y)
    {
        return a;
    }
    else if(b.y>=a.y && b.y>=c.y)
    {
        return b;
    }
    else
    {
        return c;
    }
}

homogeneous_point min_Y(homogeneous_point a, homogeneous_point b, homogeneous_point c)
{
    if(a.y<=b.y && a.y<=c.y)
    {
        return a;
    }
    else if(b.y<=a.y && b.y<=c.y)
    {
        return b;
    }
    else
    {
        return c;
    }
}

homogeneous_point mid_Y(homogeneous_point a, homogeneous_point b, homogeneous_point c)
{
    homogeneous_point temp1=max_Y(a,b,c);
    homogeneous_point temp2=min_Y(a,b,c);

    if((a==temp1 && b==temp2) || (a==temp2 && b==temp1))
    {
        return c;
    }
    else if((a==temp1 && c==temp2) || (a==temp2 && c==temp1))
    {
        return b;
    }
    else{
        return a;
    }
}



double eye_x, eye_y, eye_z;
double look_x, look_y, look_z;
double up_x, up_y, up_z;
double fov_x, fov_y, aspectRatio, near, far;
color backgroud;
int screen_x, screen_y;

int tot_tri=0;
int tot_tri1=0;
int tot_tri2=0;
color Tri_color1[1000];
color Tri_color2[1000];
color Tri_color3[1000];


void scan_convert() {
    ifstream stage3;
    stage3.open("stage3.txt");

    //color** pixels = new color*[screen_x];
    //double** zs = new double*[screen_x];

    //color pixels[screen_x+1][screen_y+1];
    //double zs[screen_x+1][screen_y+1];

    vector<vector<color>>pixels(screen_x+1,vector<color>(screen_y+1));
    vector<vector<double>>zs(screen_x+1,vector<double>(screen_y+1));

    for (int i = 0; i < screen_x; i++) {
        //pixels[i] = new color [screen_y];
        for (int j = 0; j < screen_y; j++) {
            pixels[i][j] = backgroud;
        }
        //zs[i] = new double [screen_y];
        for (int j = 0; j < screen_y; j++) {
            zs[i][j] = 20; // a very large value intended as +INFINITY
        }
    }

    cout<<tot_tri2<<endl;

    for(int ii=0;ii<tot_tri2;ii++)
    {
        //cout<<"edi kiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii "<<ii+1<<endl;
        homogeneous_point hp_arra[4];
        for(int jj=0;jj<3;jj++){
            double xx,yy,zz;
            stage3 >> xx >> yy >> zz;
            homogeneous_point hpp1(xx,yy,zz);
            hp_arra[jj]=hpp1;
            //homogeneous_point hp2=temp_mat*hp1;
            //stage3<< hp2.x <<" "<<hp2.y<<" "<< hp2.z<<endl;
        }
        homogeneous_point mxY=max_Y(hp_arra[0],hp_arra[1],hp_arra[2]);
        homogeneous_point mdY=mid_Y(hp_arra[0],hp_arra[1],hp_arra[2]);
        homogeneous_point mnY=min_Y(hp_arra[0],hp_arra[1],hp_arra[2]);

		if(mxY.y==mdY.y)
        {
            double point1=floor(((double)screen_y*(1-mnY.y)-1.0)/2.0);
            int pnt1=(int)point1;//mid/min

            double point2=ceil(((double)screen_y*(1-mxY.y)-1.0)/2.0);
            int pnt2=(int)point2;//max

            if(mxY.y<-1 || mnY.y>1)
            {
                pnt2=-1;
                pnt1=-2;

            }
            else if(mxY.y>1 && (mnY.y<1 && mnY.y>-1))
            {
                pnt2=0;
                pnt1=pnt1;
            }
            else if((mxY.y<1 && mxY.y>-1) && mnY.y<-1)
            {
                pnt2=pnt2;
                pnt1=screen_y-1;
            }
            else if(mxY.y>1 && mnY.y<-1)
            {
                pnt2=0;
                pnt1=screen_y-1;
            }
            int ii1;
            for(ii1=pnt2;ii1<=pnt1;ii1++)
            {

                double y_val=1.0-(2.0*ii1+1.0)/(double)screen_y;
                homogeneous_point hp1(0.0,0.0,0.0);
                homogeneous_point hp2(0.0,0.0,0.0);

                if(mnY.y==mdY.y || mxY.y==mnY.y)
                {
                    hp1.x=1.2;
                    hp2.x=1.1;
                }
                else{
                    hp1=Line_Plane_Intersect1(mdY,mnY,y_val);
                    hp2=Line_Plane_Intersect1(mnY,mxY,y_val);
                }

                double mxX,mnX,mxZ,mnZ;
                if(hp1.x>=hp2.x)
                {
                    mxX=hp1.x;
                    mnX=hp2.x;

                    mxZ=hp1.z;
                    mnZ=hp2.z;

                }
                else{
                    mxX=hp2.x;
                    mnX=hp1.x;

                    mxZ=hp2.z;
                    mnZ=hp1.z;
                }
                double pntX1=(mnX+1.0)*(double)screen_x/2.0;
                double pntX2=(mxX+1.0)*(double)screen_x/2.0;

                int px1=(int)pntX1;
                int px2=(int)pntX2;
                int jj1;
                for(jj1=px1;jj1<=px2;jj1++)
                {
                    if(jj1>=0 && jj1<=screen_x)
                    {
                        double xp=-1.0+(2.0*(double)jj1)/screen_x;
                        double zp=mxZ-(mxZ-mnZ)*(mxX-xp)/(mxX-mnX);
                        if(zs[jj1][ii1]>=zp+epsilon)
                        {
                            zs[jj1][ii1]=zp;
                            pixels[jj1][ii1]=Tri_color3[ii];
                        }
                    }
                }


            }


        }
        else if(mdY.y==mnY.y)
        {
            double point1=floor(((double)screen_y*(1-mdY.y)-1.0)/2.0);
            int pnt1=(int)point1;//mid/min

            double point2=ceil(((double)screen_y*(1-mxY.y)-1.0)/2.0);
            int pnt2=(int)point2;//max

            if(mxY.y<-1 || mdY.y>1)
            {
                pnt2=-1;
                pnt1=-2;

            }
            else if(mxY.y>1 && (mdY.y<1 && mdY.y>-1))
            {
                pnt2=0;
                pnt1=pnt1;
            }
            else if((mxY.y<1 && mxY.y>-1) && mdY.y<-1)
            {
                pnt2=pnt2;
                pnt1=screen_y-1;
            }
            else if(mxY.y>1 && mdY.y<-1)
            {
                pnt2=0;
                pnt1=screen_y-1;
            }
            int ii1;
            for(ii1=pnt2;ii1<=pnt1;ii1++)
            {

                double y_val=1.0-(2.0*ii1+1.0)/(double)screen_y;

                homogeneous_point hp1(0.0,0.0,0.0);
                homogeneous_point hp2(0.0,0.0,0.0);

                if(mdY.y==mxY.y || mnY.y==mxY.y)
                {
                    hp1.x=1.2;
                    hp2.x=1.1;
                }
                else
                {
                    hp1=Line_Plane_Intersect1(mdY,mxY,y_val);
                    hp2=Line_Plane_Intersect1(mnY,mxY,y_val);
                }

                double mxX,mnX,mxZ,mnZ;
                if(hp1.x>=hp2.x)
                {
                    mxX=hp1.x;
                    mnX=hp2.x;

                    mxZ=hp1.z;
                    mnZ=hp2.z;

                }
                else{
                    mxX=hp2.x;
                    mnX=hp1.x;

                    mxZ=hp2.z;
                    mnZ=hp1.z;
                }
                double pntX1=(mnX+1.0)*(double)screen_x/2.0;
                double pntX2=(mxX+1.0)*(double)screen_x/2.0;

                int px1=(int)pntX1;
                int px2=(int)pntX2;
                int jj1;
                for(jj1=px1;jj1<=px2;jj1++)
                {
                    if(jj1>=0 && jj1<=screen_x)
                    {
                        double xp=-1.0+(2.0*(double)jj1)/screen_x;
                        double zp=mxZ-(mxZ-mnZ)*(mxX-xp)/(mxX-mnX);
                        if(zs[jj1][ii1]>=zp+epsilon)
                        {
                            zs[jj1][ii1]=zp;
                            pixels[jj1][ii1]=Tri_color3[ii];
                        }
                    }
                }


            }


        }
        else
        {
            double point1=floor(((double)screen_y*(1-mnY.y)-1.0)/2.0);
            int pnt1=(int)point1;//mid/min

            double point2=ceil(((double)screen_y*(1-mxY.y)-1.0)/2.0);
            int pnt2=(int)point2;//max

            int ii1;
            for(ii1=pnt2;ii1<=pnt1;ii1++)
            {
                if(ii1>=0 && ii1<=screen_y){

                    double y_val=1.0-(2.0*ii1+1.0)/(double)screen_y;
                    homogeneous_point hp1(0.0,0.0,0.0);
                    homogeneous_point hp2(0.0,0.0,0.0);

                    if(y_val>=mdY.y){
						if(mxY.y==mdY.y || mxY.y==mnY.y){
							hp1.x=1.2;
							hp2.x=1.1;
						}
						else
						{
							hp1=Line_Plane_Intersect1(mdY,mxY,y_val);
							hp2=Line_Plane_Intersect1(mnY,mxY,y_val);
						}
					}
					else
					{
						if(mnY.y==mdY.y || mxY.y==mnY.y)
						{
							hp1.x=1.2;
							hp2.x=1.1;
						}
						else{
							hp1=Line_Plane_Intersect1(mdY,mnY,y_val);
							hp2=Line_Plane_Intersect1(mnY,mxY,y_val);

						}

					}

                    double mxX,mnX,mxZ,mnZ;
                    if(hp1.x>=hp2.x)
                    {
                        mxX=hp1.x;
                        mnX=hp2.x;

                        mxZ=hp1.z;
                        mnZ=hp2.z;

                    }
                    else{
                        mxX=hp2.x;
                        mnX=hp1.x;

                        mxZ=hp2.z;
                        mnZ=hp1.z;
                    }
                    double pntX1=(mnX+1.0)*(double)screen_x/2.0;
                    double pntX2=(mxX+1.0)*(double)screen_x/2.0;

                    int px1=(int)pntX1;
                    int px2=(int)pntX2;
                    int jj1;
                    for(jj1=px1;jj1<=px2;jj1++)
                    {
                        if(jj1>=0 && jj1<=screen_x)
                        {
                            double xp=-1.0+(2.0*(double)jj1)/screen_x;
                            double zp=mxZ-(mxZ-mnZ)*(mxX-xp)/(mxX-mnX);
                            if(zs[jj1][ii1]>=zp+epsilon)
                            {
                                zs[jj1][ii1]=zp;
                                pixels[jj1][ii1]=Tri_color3[ii];
                            }
                        }
                    }


                }


            }

        }
    }

    // perform scan conversion, populate the 2D array pixels
    // the array zs is the z-buffer.


    // the following code generates a bmp image. do not change this.
    bitmap_image image(screen_x, screen_y);
    for (int x = 0; x < screen_x; x++) {
        for (int y = 0; y < screen_y; y++) {
            image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
        }
    }


    cout<<"Ashe"<<endl;
    image.save_image("out.bmp");

    // free the dynamically allocated memory

}


void stage3()
{
    if (near == far) return;
    ifstream stage2;
    ofstream temp1;
    ifstream temp2;
    ofstream stage3;
    stage2.open ("stage2.txt");
    temp1.open("temp1.txt");
    stage3.open ("stage3.txt");
    stage3 << std::fixed;
    stage3 << std::setprecision(7);

    temp1 << std::fixed;
    temp1 << std::setprecision(7);



    fov_x=fov_y*aspectRatio;
    double t= near * tan(fov_y*2*pi/(360*2));
    double r= near * tan(fov_x*2*pi/(360*2));

    matrix temp_mat(4);
    temp_mat=temp_mat.Zero();

    temp_mat.values[0][0]=near/r;
    temp_mat.values[1][1]=near/t;
    temp_mat.values[2][2]=-((far+near)/(far-near));
    temp_mat.values[3][2]=-1;
    temp_mat.values[2][3]=-((2*far*near)/(far-near));

    // process input from stage2 and write to stage3
    for(int ii=0;ii<tot_tri;ii++)
    {
        homogeneous_point hp_arra[4];
        for(int jj=0;jj<3;jj++){
            double xx,yy,zz;
            stage2 >> xx >> yy >> zz;
            homogeneous_point hp1(xx,yy,zz);
            hp_arra[jj]=hp1;
            //homogeneous_point hp2=temp_mat*hp1;
            //stage3<< hp2.x <<" "<<hp2.y<<" "<< hp2.z<<endl;
        }

        //stage3<<endl;
        //cout<< xx << " "<<yy<<" "<<zz<<endl;
        double zz0=hp_arra[0].z;
        double zz1=hp_arra[1].z;
        double zz2=hp_arra[2].z;

        if(hp_arra[0].z<-far && hp_arra[1].z<-far && hp_arra[2].z<-far)
        {

        }
        else if(hp_arra[0].z>=-far && hp_arra[1].z>=-far && hp_arra[2].z>=-far)
        {
            temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            temp1<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];
        }
        //touch condition
        else if(hp_arra[0].z==-far)
        {
            if((hp_arra[1].z<-far && hp_arra[2].z>-far))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[1],hp_arra[2],-far);

                temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                temp1<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                temp1<<endl;

                tot_tri1++;
                Tri_color2[tot_tri1-1]=Tri_color1[ii];

            }
            else if((hp_arra[1].z>-far && hp_arra[2].z<-far))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[1],hp_arra[2],-far);

                temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                temp1<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                temp1<<endl;

                tot_tri1++;
                Tri_color2[tot_tri1-1]=Tri_color1[ii];
            }
            else if(hp_arra[1].z==-far)
            {
                if(hp_arra[2].z>-far)
                {
                    temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    temp1<<endl;

                    tot_tri1++;
                    Tri_color2[tot_tri1-1]=Tri_color1[ii];
                }
                else if(hp_arra[2].z<-far)
                {

                }
            }
            else if(hp_arra[2].z==-far)
            {
                if(hp_arra[1].z>-far)
                {
                    temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    temp1<<endl;

                    tot_tri1++;
                    Tri_color2[tot_tri1-1]=Tri_color1[ii];
                }
                else if(hp_arra[1].z<-far)
                {

                }
            }
            else if((hp_arra[1].z>-far && hp_arra[2].z>-far))
            {

                temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                temp1<<endl;

                tot_tri1++;
                Tri_color2[tot_tri1-1]=Tri_color1[ii];

            }
            else if((hp_arra[1].z<-far && hp_arra[2].z<-far))
            {

            }
        }

        else if(hp_arra[1].z==-far)
        {
            if((hp_arra[0].z<-far && hp_arra[2].z>-far))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[0],hp_arra[2],-far);

                temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                temp1<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                temp1<<endl;

                tot_tri1++;
                Tri_color2[tot_tri1-1]=Tri_color1[ii];

            }
            else if((hp_arra[0].z>-far && hp_arra[2].z<-far))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[0],hp_arra[2],-far);

                temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                temp1<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                temp1<<endl;

                tot_tri1++;
                Tri_color2[tot_tri1-1]=Tri_color1[ii];
            }
            else if(hp_arra[0].z==-far)
            {
                if(hp_arra[2].z>-far)
                {
                    temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    temp1<<endl;

                    tot_tri1++;
                    Tri_color2[tot_tri1-1]=Tri_color1[ii];
                }
                else if(hp_arra[2].z<-far)
                {

                }
            }
            else if(hp_arra[2].z==-far)
            {
                if(hp_arra[0].z>-far)
                {
                    temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    temp1<<endl;

                    tot_tri1++;
                    Tri_color2[tot_tri1-1]=Tri_color1[ii];
                }
                else if(hp_arra[0].z<-far)
                {

                }
            }
            else if((hp_arra[0].z>-far && hp_arra[2].z>-far))
            {

                temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                temp1<<endl;

                tot_tri1++;
                Tri_color2[tot_tri1-1]=Tri_color1[ii];

            }
            else if((hp_arra[0].z<-far && hp_arra[2].z<-far))
            {

            }
        }

        else if(hp_arra[2].z==-far)
        {
            if((hp_arra[1].z<-far && hp_arra[0].z>-far))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[1],hp_arra[0],-far);

                temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                temp1<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                temp1<<endl;

                tot_tri1++;
                Tri_color2[tot_tri1-1]=Tri_color1[ii];

            }
            else if((hp_arra[1].z>-far && hp_arra[0].z<-far))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[1],hp_arra[0],-far);

                temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                temp1<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                temp1<<endl;

                tot_tri1++;
                Tri_color2[tot_tri1-1]=Tri_color1[ii];
            }
            else if(hp_arra[1].z==-far)
            {
                if(hp_arra[0].z>-far)
                {
                    temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    temp1<<endl;

                    tot_tri1++;
                    Tri_color2[tot_tri1-1]=Tri_color1[ii];
                }
                else if(hp_arra[0].z<-far)
                {

                }
            }
            else if(hp_arra[0].z==-far)
            {
                if(hp_arra[1].z>-far)
                {
                    temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    temp1<<endl;

                    tot_tri1++;
                    Tri_color2[tot_tri1-1]=Tri_color1[ii];
                }
                else if(hp_arra[1].z<-far)
                {

                }
            }
            else if((hp_arra[1].z>-far && hp_arra[0].z>-far))
            {

                temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                temp1<<endl;

                tot_tri1++;
                Tri_color2[tot_tri1-1]=Tri_color1[ii];

            }
            else if((hp_arra[1].z<-far && hp_arra[0].z<-far))
            {

            }
        }
        //ekta vitre, 2 ta baire
        else if(hp_arra[0].z>-far && hp_arra[1].z<-far && hp_arra[2].z<-far)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[0],hp_arra[1],-far);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[0],hp_arra[2],-far);
            temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            temp1<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            temp1<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            temp1<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];


        }
        else if(hp_arra[1].z>-far && hp_arra[2].z<-far && hp_arra[0].z<-far)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[1],hp_arra[2],-far);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[1],hp_arra[0],-far);
            temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            temp1<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            temp1<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            temp1<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];

        }
        else if(hp_arra[2].z>-far && hp_arra[0].z<-far && hp_arra[1].z<-far)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[2],hp_arra[1],-far);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[2],hp_arra[0],-far);
            temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            temp1<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            temp1<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            temp1<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];

        }
        //2 ta vitre, ekta baire
        else if(hp_arra[0].z>-far && hp_arra[1].z>-far && hp_arra[2].z<-far)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[0],hp_arra[2],-far);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[1],hp_arra[2],-far);
            temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            temp1<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];

            temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            temp1<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            temp1<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            temp1<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];


        }
        else if(hp_arra[1].z>-far && hp_arra[2].z>-far && hp_arra[0].z<-far)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[1],hp_arra[0],-far);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[2],hp_arra[0],-far);
            temp1<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            temp1<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];

            temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            temp1<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            temp1<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            temp1<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];

        }
        else if(hp_arra[2].z>-far && hp_arra[0].z>-far && hp_arra[1].z<-far)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[2],hp_arra[1],-far);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[0],hp_arra[1],-far);
            temp1<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            temp1<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];


            temp1<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            temp1<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            temp1<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            temp1<<endl;

            tot_tri1++;
            Tri_color2[tot_tri1-1]=Tri_color1[ii];

        }
    }
    temp1.close();


    temp2.open("temp1.txt");
    for(int ii=0;ii<tot_tri1;ii++)
    {
        homogeneous_point hp_arra[4];
        for(int jj=0;jj<3;jj++){
            double xx,yy,zz;
            temp2 >> xx >> yy >> zz;
            homogeneous_point hp1(xx,yy,zz);
            hp_arra[jj]=hp1;
            //homogeneous_point hp2=temp_mat*hp1;
            //stage3<< hp2.x <<" "<<hp2.y<<" "<< hp2.z<<endl;
        }
        //stage3<<endl;
        //cout<< xx << " "<<yy<<" "<<zz<<endl;
        double zz0=hp_arra[0].z;
        double zz1=hp_arra[1].z;
        double zz2=hp_arra[2].z;

        if(hp_arra[0].z>-near && hp_arra[1].z>-near && hp_arra[2].z>-near)
        {

        }
        else if(hp_arra[0].z<=-near && hp_arra[1].z<=-near && hp_arra[2].z<=-near)
        {
            hp_arra[0]=temp_mat*hp_arra[0];
            hp_arra[1]=temp_mat*hp_arra[1];
            hp_arra[2]=temp_mat*hp_arra[2];

            stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            stage3<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];
        }
        //touch condition
        else if(hp_arra[0].z==-near)
        {
            if((hp_arra[1].z>-near && hp_arra[2].z<-near))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[1],hp_arra[2],-near);

                hp_arra[0]=temp_mat*hp_arra[0];
                tmpp=temp_mat*tmpp;
                hp_arra[2]=temp_mat*hp_arra[2];

                stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                stage3<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                stage3<<endl;

                tot_tri2++;
                Tri_color3[tot_tri2-1]=Tri_color2[ii];

            }
            else if((hp_arra[1].z<-near && hp_arra[2].z>-near))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[1],hp_arra[2],-near);

                hp_arra[0]=temp_mat*hp_arra[0];
                hp_arra[1]=temp_mat*hp_arra[1];
                tmpp=temp_mat*tmpp;

                stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                stage3<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                stage3<<endl;

                tot_tri2++;
                Tri_color3[tot_tri2-1]=Tri_color2[ii];
            }
            else if(hp_arra[1].z==-near)
            {
                if(hp_arra[2].z<-near)
                {
                    hp_arra[0]=temp_mat*hp_arra[0];
                    hp_arra[1]=temp_mat*hp_arra[1];
                    hp_arra[2]=temp_mat*hp_arra[2];

                    stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    stage3<<endl;

                    tot_tri2++;
                    Tri_color3[tot_tri2-1]=Tri_color2[ii];
                }
                else if(hp_arra[2].z>-near)
                {

                }
            }
            else if(hp_arra[2].z==-near)
            {
                if(hp_arra[1].z<-near)
                {
                    hp_arra[0]=temp_mat*hp_arra[0];
                    hp_arra[1]=temp_mat*hp_arra[1];
                    hp_arra[2]=temp_mat*hp_arra[2];

                    stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    stage3<<endl;

                    tot_tri2++;
                    Tri_color3[tot_tri2-1]=Tri_color2[ii];
                }
                else if(hp_arra[1].z>-near)
                {

                }
            }
            else if((hp_arra[1].z<-near && hp_arra[2].z<-near))
            {

                hp_arra[0]=temp_mat*hp_arra[0];
                hp_arra[1]=temp_mat*hp_arra[1];
                hp_arra[2]=temp_mat*hp_arra[2];

                stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                stage3<<endl;

                tot_tri2++;
                Tri_color3[tot_tri2-1]=Tri_color2[ii];

            }
            else if((hp_arra[1].z>-near && hp_arra[2].z>-near))
            {

            }
        }

        else if(hp_arra[1].z==-near)
        {
            if((hp_arra[0].z>-near && hp_arra[2].z<-near))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[0],hp_arra[2],-near);

                hp_arra[1]=temp_mat*hp_arra[1];
                tmpp=temp_mat*tmpp;
                hp_arra[2]=temp_mat*hp_arra[2];

                stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                stage3<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                stage3<<endl;

                tot_tri2++;
                Tri_color3[tot_tri2-1]=Tri_color2[ii];

            }
            else if((hp_arra[0].z<-near && hp_arra[2].z>-near))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[0],hp_arra[2],-near);

                hp_arra[1]=temp_mat*hp_arra[1];
                tmpp=temp_mat*tmpp;
                hp_arra[0]=temp_mat*hp_arra[0];

                stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                stage3<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                stage3<<endl;

                tot_tri2++;
                Tri_color3[tot_tri2-1]=Tri_color2[ii];
            }
            else if(hp_arra[0].z==-near)
            {
                if(hp_arra[2].z<-near)
                {
                    hp_arra[0]=temp_mat*hp_arra[0];
                    hp_arra[1]=temp_mat*hp_arra[1];
                    hp_arra[2]=temp_mat*hp_arra[2];

                    stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    stage3<<endl;

                    tot_tri2++;
                    Tri_color3[tot_tri2-1]=Tri_color2[ii];
                }
                else if(hp_arra[2].z>-near)
                {

                }
            }
            else if(hp_arra[2].z==-near)
            {
                if(hp_arra[0].z<-near)
                {
                    hp_arra[0]=temp_mat*hp_arra[0];
                    hp_arra[1]=temp_mat*hp_arra[1];
                    hp_arra[2]=temp_mat*hp_arra[2];

                    stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    stage3<<endl;

                    tot_tri2++;
                    Tri_color3[tot_tri2-1]=Tri_color2[ii];
                }
                else if(hp_arra[0].z>-near)
                {

                }
            }
            else if((hp_arra[0].z<-near && hp_arra[2].z<-near))
            {
                hp_arra[0]=temp_mat*hp_arra[0];
                hp_arra[1]=temp_mat*hp_arra[1];
                hp_arra[2]=temp_mat*hp_arra[2];

                stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                stage3<<endl;

                tot_tri2++;
                Tri_color3[tot_tri2-1]=Tri_color2[ii];

            }
            else if((hp_arra[0].z>-near && hp_arra[2].z>-near))
            {

            }
        }

        else if(hp_arra[2].z==-near)
        {
            if((hp_arra[1].z>-near && hp_arra[0].z<-near))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[1],hp_arra[0],-near);

                hp_arra[0]=temp_mat*hp_arra[0];
                tmpp=temp_mat*tmpp;
                hp_arra[2]=temp_mat*hp_arra[2];

                stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                stage3<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                stage3<<endl;

                tot_tri2++;
                Tri_color3[tot_tri2-1]=Tri_color2[ii];

            }
            else if((hp_arra[1].z<-near && hp_arra[0].z>-near))
            {
                homogeneous_point tmpp=Line_Plane_Intersect(hp_arra[1],hp_arra[0],-near);

                tmpp=temp_mat*tmpp;
                hp_arra[1]=temp_mat*hp_arra[1];
                hp_arra[2]=temp_mat*hp_arra[2];

                stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                stage3<< tmpp.x<<" "<<tmpp.y<<" "<<tmpp.z<<endl;
                stage3<<endl;

                tot_tri2++;
                Tri_color3[tot_tri2-1]=Tri_color2[ii];
            }
            else if(hp_arra[1].z==-near)
            {
                if(hp_arra[0].z<-near)
                {
                    hp_arra[0]=temp_mat*hp_arra[0];
                    hp_arra[1]=temp_mat*hp_arra[1];
                    hp_arra[2]=temp_mat*hp_arra[2];

                    stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    stage3<<endl;

                    tot_tri2++;
                    Tri_color3[tot_tri2-1]=Tri_color2[ii];
                }
                else if(hp_arra[0].z>-near)
                {

                }
            }
            else if(hp_arra[0].z==-near)
            {
                if(hp_arra[1].z<-near)
                {
                    hp_arra[0]=temp_mat*hp_arra[0];
                    hp_arra[1]=temp_mat*hp_arra[1];
                    hp_arra[2]=temp_mat*hp_arra[2];

                    stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                    stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                    stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                    stage3<<endl;

                    tot_tri2++;
                    Tri_color3[tot_tri2-1]=Tri_color2[ii];
                }
                else if(hp_arra[1].z>-near)
                {

                }
            }
            else if((hp_arra[1].z<-near && hp_arra[0].z<-near))
            {
                hp_arra[0]=temp_mat*hp_arra[0];
                hp_arra[1]=temp_mat*hp_arra[1];
                hp_arra[2]=temp_mat*hp_arra[2];

                stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
                stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
                stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
                stage3<<endl;

                tot_tri2++;
                Tri_color3[tot_tri2-1]=Tri_color2[ii];

            }
            else if((hp_arra[1].z>-near && hp_arra[0].z>-near))
            {

            }
        }
        //ekta vitre, 2 ta baire
        else if(hp_arra[0].z<-near && hp_arra[1].z>-near && hp_arra[2].z>-near)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[0],hp_arra[1],-near);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[0],hp_arra[2],-near);

            hp_arra[0]=temp_mat*hp_arra[0];
            hp1=temp_mat*hp1;
            hp2=temp_mat*hp2;

            stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            stage3<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            stage3<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            stage3<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];


        }
        else if(hp_arra[1].z<-near && hp_arra[2].z>-near && hp_arra[0].z>-near)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[1],hp_arra[2],-near);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[1],hp_arra[0],-near);

            hp_arra[1]=temp_mat*hp_arra[1];
            hp1=temp_mat*hp1;
            hp2=temp_mat*hp2;


            stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            stage3<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            stage3<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            stage3<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];

        }
        else if(hp_arra[2].z<-near && hp_arra[0].z>-near && hp_arra[1].z>-near)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[2],hp_arra[1],-near);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[2],hp_arra[0],-near);

            hp_arra[2]=temp_mat*hp_arra[2];
            hp1=temp_mat*hp1;
            hp2=temp_mat*hp2;

            stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            stage3<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            stage3<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            stage3<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];

        }
        //2 ta vitre, ekta baire
        else if(hp_arra[0].z<-near && hp_arra[1].z<-near && hp_arra[2].z>-near)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[0],hp_arra[2],-near);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[1],hp_arra[2],-near);

            hp_arra[0]=temp_mat*hp_arra[0];
            hp_arra[1]=temp_mat*hp_arra[1];
            hp1=temp_mat*hp1;
            hp2=temp_mat*hp2;

            stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            stage3<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];

            stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            stage3<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            stage3<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            stage3<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];


        }
        else if(hp_arra[1].z<-near && hp_arra[2].z<-near && hp_arra[0].z>-near)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[1],hp_arra[0],-near);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[2],hp_arra[0],-near);

            hp_arra[2]=temp_mat*hp_arra[2];
            hp_arra[1]=temp_mat*hp_arra[1];
            hp1=temp_mat*hp1;
            hp2=temp_mat*hp2;

            stage3<< hp_arra[1].x<<" "<<hp_arra[1].y<<" "<<hp_arra[1].z<<endl;
            stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            stage3<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];

            stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            stage3<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            stage3<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            stage3<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];

        }
        else if(hp_arra[2].z<-near && hp_arra[0].z<-near && hp_arra[1].z>-near)
        {
            homogeneous_point hp1=Line_Plane_Intersect(hp_arra[2],hp_arra[1],-near);
            homogeneous_point hp2=Line_Plane_Intersect(hp_arra[0],hp_arra[1],-near);

            hp_arra[0]=temp_mat*hp_arra[0];
            hp_arra[2]=temp_mat*hp_arra[2];
            hp1=temp_mat*hp1;
            hp2=temp_mat*hp2;

            stage3<< hp_arra[2].x<<" "<<hp_arra[2].y<<" "<<hp_arra[2].z<<endl;
            stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            stage3<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];

            stage3<< hp_arra[0].x<<" "<<hp_arra[0].y<<" "<<hp_arra[0].z<<endl;
            stage3<< hp1.x<<" "<<hp1.y<<" "<<hp1.z<<endl;
            stage3<< hp2.x<<" "<<hp2.y<<" "<<hp2.z<<endl;
            stage3<<endl;

            tot_tri2++;
            Tri_color3[tot_tri2-1]=Tri_color2[ii];

        }
    }
    temp2.close();


    stage3.close();
    stage2.close();

}

void stage2()
{
    ifstream stage1;
    ofstream stage2;
    stage1.open ("stage1.txt");
    stage2.open ("stage2.txt");
    stage2 << std::fixed;
    stage2 << std::setprecision(7);

    // collect input from stage1 and process, write output to stage2
    Vector eye(eye_x,eye_y,eye_z);
    Vector look(look_x,look_y,look_z);
    Vector up(up_x,up_y,up_z);

    Vector l = look - eye;
    l.normalize();

    Vector r(0.0,0.0,0.0);
    r=r.cross(l,up);
    r.normalize();

    Vector u(0.0,0.0,0.0);
    u=u.cross(r,l);

    matrix T(4);
    T=T.make_identity(4);

    T.values[0][3]=-eye_x;
    T.values[1][3]=-eye_y;
    T.values[2][3]=-eye_z;

    matrix R1(4);
    R1=R1.make_identity(4);

    R1.values[0][0]=r.x;
    R1.values[0][1]=r.y;
    R1.values[0][2]=r.z;

    R1.values[1][0]=u.x;
    R1.values[1][1]=u.y;
    R1.values[1][2]=u.z;

    R1.values[2][0]=-l.x;
    R1.values[2][1]=-l.y;
    R1.values[2][2]=-l.z;

    matrix V=R1*T;
    //V.print();

    for(int ii=0;ii<tot_tri;ii++)
    {
        for(int jj=0;jj<3;jj++){
            double xx,yy,zz;
            stage1 >> xx >> yy >> zz;
            homogeneous_point hp1(xx,yy,zz);
            homogeneous_point hp2=V*hp1;
            stage2<< hp2.x <<" "<<hp2.y<<" "<< hp2.z<<endl;
        }
        stage2<<endl;
        //cout<< xx << " "<<yy<<" "<<zz<<endl;
    }


    stage1.close();
    stage2.close();

}

void stage1()
{
    ifstream scene;
    ofstream stage1;
    scene.open ("scene.txt");
    stage1.open ("stage1.txt");
    stage1 << std::fixed;
    stage1 << std::setprecision(7);

    string command;

    scene >> eye_x >> eye_y >> eye_z;
    scene >> look_x >> look_y >> look_z;
    scene >> up_x >> up_y >> up_z;
    scene >> fov_y >> aspectRatio >> near >> far;
    scene >> screen_x >> screen_y;
    scene >> backgroud.r >> backgroud.g >> backgroud.b;

    // take other commands as input from scene in a loop
    std::stack<matrix> stk;
    matrix mat1(4);
    mat1=mat1.make_identity(4);
    //mat1.print();
    stk.push(mat1);

    int cnt=0;
    int stk_size_bfr_push;
    int arra[1000];
    memset(arra,0,sizeof(arra));


    while(1==1)
    {
        scene >> command;
        if(command=="triangle")
        {
            tot_tri++;
            for(int ii=1;ii<=3;ii++)
            {
                double triX,triY,triZ;
                scene >> triX >> triY >> triZ;
                homogeneous_point hp1(triX, triY, triZ);
                homogeneous_point hp2=stk.top()*hp1;
                stage1<< hp2.x <<" "<<hp2.y<<" "<< hp2.z<<endl;

            }
            stage1<<endl;
            int colR, colG, colB;
            scene >> colR >> colG >> colB;
            color cl(colR,colG,colB);
            Tri_color1[tot_tri-1]=cl;
        }
        else if(command=="translate")
        {
            double tx, ty, tz;
            scene >> tx >> ty >> tz;
            matrix temp_Mat(4);
            temp_Mat=temp_Mat.translate(tx,ty,tz);
            temp_Mat=stk.top()*temp_Mat;
            stk.push(temp_Mat);
            arra[cnt]++;

        }
        else if(command=="scale")
        {
            double sx, sy, sz;
            scene >> sx >> sy >> sz;
            matrix temp_Mat(4);
            temp_Mat=temp_Mat.scale(sx,sy,sz);
            temp_Mat=stk.top()*temp_Mat;
            stk.push(temp_Mat);
            arra[cnt]++;

        }
        else if(command=="rotate")
        {
            double angle1, ax, ay, az;
            scene >> angle1 >> ax >> ay >> az;

            Vector a_of_Rot(ax,ay,az);
            a_of_Rot.normalize();

            Vector c1(1.0,0.0,0.0);//i
            Vector c2(0.0,1.0,0.0);//j
            Vector c3(0.0,0.0,1.0);//k

            c1=c1.R(c1,a_of_Rot,angle1);
            c2=c2.R(c2,a_of_Rot,angle1);
            c3=c3.R(c3,a_of_Rot,angle1);

            matrix temp_Mat(4);
            temp_Mat=temp_Mat.make_identity(4);

            temp_Mat.values[0][0]=c1.x;
            temp_Mat.values[1][0]=c1.y;
            temp_Mat.values[2][0]=c1.z;

            temp_Mat.values[0][1]=c2.x;
            temp_Mat.values[1][1]=c2.y;
            temp_Mat.values[2][1]=c2.z;

            temp_Mat.values[0][2]=c3.x;
            temp_Mat.values[1][2]=c3.y;
            temp_Mat.values[2][2]=c3.z;

            temp_Mat=stk.top()*temp_Mat;
            stk.push(temp_Mat);
            arra[cnt]++;


        }
        else if(command=="push")
        {
            cnt++;

        }
        else if(command=="pop")
        {
            if(cnt!=0){
                while(arra[cnt]>0)
                {
                    //cout<<"DHukse"<<endl;
                    stk.pop();
                    arra[cnt]--;
                }
                arra[cnt]=0;
                cnt--;
            }

        }
        else if(command=="end")
        {
            break;
        }
    }
    // process accordingly
    // write to stage1

    scene.close();
    stage1.close();

}

int main()
{
    cout << std::fixed;
    cout << std::setprecision(4);

    stage1();
    stage2();
    stage3();
    scan_convert();

    return 0;
}
