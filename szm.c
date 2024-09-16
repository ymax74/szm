#include <stdlib.h>
#include "stdio.h"
#include <math.h>
#define pi        3.1415926535897932

class cmnkB  
{
public:
	cmnkB();
	virtual ~cmnkB();
	void initmnk(int km, int nvar);
	void setel(int km, int nvar,long double el);
	void  slvls2(int s);
	long double * sgmnk(int km, int sn);
	long double* ptlst(int km, int sn);

	long double ** a;
	long double ** c;
	long double * xf;
	long double * sgm;
	int mkm,mnvar;
	double* a25,*xf25,*cu;
	double** aa;
	double* p25();



};

class IEx
{
public:
	IEx();
	~IEx();
	void addob(int ism,int nx,int ny,int xx,int xy,int Imax,int npix);
	int nbufoxy,obty;
	int* aism,*anx,*any,*axx,*axy,*aImax,*anpix;
	int* anob;
	int szbuf;
	void sortob();
	void GetCen2(unsigned short** ap17img,unsigned short ** nesosdQ,int tt100, int w, int h,int t765, int t65);
	cmnkB omnk;
	void drawLine(int x1, int y1, int x2, int y2,int*ax,int*ay,int* n);
	float *ax,*ay,*az,*zmax,*lnsh,*tetha;
	float *npx;
	int nstars;
	

	void GetSources(unsigned short* ccdimgC, int w, int h, double dS, int nLimSh, int urnas,
		int* nstr,
		float** axS,
		float** ayS,
		float** azS,
		float** zmaxS,
		float** lnshS,
		float** tethaS,
		float** npxS);
	unsigned short** nesosd;
	unsigned short** ap17img;
	long** img;
	int svh,svw,svt3;
	char asrv[2][2][2][2][2][2];
	char asrv3[2][2][2][2][2][2];
	char asrv2[2][2][2][2][2][2][2][2][2][2];
	int t3;//диаметр фильтра
	void SetDiaFil(int df);
	int xp1,yp1,xp2,yp2,xp3,yp3,xp4,yp4,xp5,yp5;
	long** GetSubBG();
	
	
};

cmnkB::cmnkB()
{
	a=NULL;
	mkm=0;
	mnvar=0;
	xf=NULL;
	sgm=NULL;
	c=NULL;
	a25=new double[25];
	xf25=new double[6];
	cu=new double[6];
	aa=new double* [25];
	int km=0;
	for (int i=0;i<25;i++) aa[i]=new double[6];
			for (int y3=-2;y3<3;y3++)
			{
				for (int x3=-2;x3<3;x3++)
				{
					aa[km][0]=x3*x3;
					aa[km][1]=x3*y3;
					aa[km][2]=y3*y3;
					aa[km][3]=x3;
					aa[km][4]=y3;
					aa[km][5]=1;
					km++;
				}
			
			}
		/*c3[0][0]=0.014285714285714285;
		c3[0][1]=-0.00000000000000000;	
		c3[0][2]=0.00000000000000000;
		c3[0][3]=-0.00000000000000000;
		c3[0][4]=0.00000000000000000;
		c3[0][5]=-0.028571428571428571;*/
	/*c4[0][0]= 0.014285714285714285;
	c4[0][5]=-0.028571428571428571;
	c4[1][1]= 0.010000000000000000;
	c4[2][2]= 0.014285714285714285;
	c4[2][5]=-0.028571428571428571;
	c4[3][3]= 0.020000000000000000;
	c4[4][4]= 0.020000000000000000;
	c4[5][0]=-0.028571428571428571;
	c4[5][2]=-0.028571428571428571;
	c4[5][5]= 0.154285714285714280;*/



}

cmnkB::~cmnkB()
{
	delete a25;
	delete xf25;
	delete cu;
	
	for (int i=0;i<25;i++) delete aa[i];
	delete aa;

	if (c!=NULL)
	{
		for(int i=0; i<mnvar+3;i++)
		{
			delete c[i];
		}
		delete c;
		c=NULL;
	}

	if (a!=NULL)
	{
		for(int i=0; i<mkm;i++)
		{
			delete a[i];
		}
		delete a;
		a=NULL;
	}

	if (xf!=NULL)
	{
		delete xf;
		xf=NULL;
	}
	
	if (sgm!=NULL)
	{
		delete sgm;
		sgm=NULL;
	}
}
void cmnkB::setel(int km, int nvar,long double el)
{
	a[km][nvar]=el;

}
void cmnkB::initmnk(int km, int nvar)
{
	if (c!=NULL)
	{
		for(int i=0; i<mnvar+3;i++)
		{
			delete c[i];
		}
		delete c;
		c=NULL;
	}

	if (a!=NULL)
	{
		for(int i=0; i<mkm;i++)
		{
			delete a[i];
		}
		delete a;
		a=NULL;
	}

	if (xf!=NULL)
	{
		delete xf;
		xf=NULL;
	}
	
	if (sgm!=NULL)
	{
		delete sgm;
		sgm=NULL;
	}
	a=new long double * [km];

	for(int i=0; i<km; i++)
	{
		a[i]=new long double [nvar+3];
	}

	c=new  long double * [nvar+3];
	int i=0;
	for( i=0; i<nvar+3; i++)
	{
		c[i]=new long double [nvar+3];
	}

	xf=new long double [nvar+3];
	sgm= new long double [nvar+3];
	mnvar=nvar;
	mkm=km;
}



void  cmnkB::slvls2(int s)
{
	
   int i,k,m;
   //a:matrls;

      for(m=0;m<s-1; m++)
      {
            for(i=m+1;i<s;i++) 
            {
                 for (k=m+1; k<s+1; k++)
                 {
                       c[i][k]=c[i][k]-c[i][m]*c[m][k]/c[m][m];
                 }
            }
      }
      for (m=s-1;m>-1; m--)
      {
           xf[m]=c[m][s];
           for (k=m+1;k<s; k++){ xf[m]=xf[m]-xf[k]*c[m][k];}
           xf[m]=xf[m]/c[m][m];

      }


}
long double *  cmnkB::sgmnk(int km,int sn)
{
	//CFile

	long double sum,sum2;
    int i,j;
     sum=0;
     for(i=0;i<km;i++)
     {
          sum2=0;
          for (j=0;j<sn; j++)
		  {
               sum2=sum2+a[i][j]*xf[j];
          };
          sum=sum+(a[i][sn]-sum2)*(a[i][sn]-sum2);
	 };
     //{d:=det(c,sn);}
     for (j=0;j<sn; j++)
     {
          sgm[j]=sqrt(sum/((km-sn)*c[j][j]));
      
     };
     sgm[sn]=sqrt(sum/((km-sn)));
	 return sgm;
}

double** Transpone(double** mas, int rows, int cols)
{
  double** rez;
  rez = new double*[cols];
  for (int i = 0; i < cols; i++)
  {
    rez[i] = new double [rows];
    for (int j = 0; j < rows; j++)
      rez[i][j] = mas[j][i];
  }
  return rez;
}


void Free(double** mas, int rows)
{
  if (mas == 0) return;    
  for (int i = 0; i < rows; i++)
    delete(mas[i]);
  delete (mas);
}


double** GetMatr(double** mas, int rows, int cols, int row, int col) {
  int di, dj;
  double** p = new double*[rows - 1];
  di = 0;
  for (int i = 0; i < rows - 1; i++) { 
    if (i == row)  
      di = 1;    
    dj = 0;
    p[i] = new double[cols - 1];
    for (int j = 0; j < cols - 1; j++) { 
      if (j == col)  
        dj = 1;    
      p[i][j] = mas[i + di][j + dj];
    }
  }
  return p;
}


double Determinant(double** mas, int m) {
  int k;
  double** p = 0;
  double d = 0;
  k = 1; 
  if (m < 1) { 
  return 0; }
  if (m == 1) { d = mas[0][0]; return(d); }
  if (m == 2) { d = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]); return(d); }
  if (m > 2) {
    for (int i = 0; i < m; i++) {
      p = GetMatr(mas, m, m, i, 0);
      d = d + k * mas[i][0] * Determinant(p, m - 1);
      k = -k;
    }
  }
  Free(p, m-1);
  return(d);
}


double** Mreverse(double** mas, int m)
{
  double** rez = new double*[m];
  double det = Determinant(mas, m); 
  for (int i = 0; i < m; i++)
  {
    rez[i] = new double[m];
    for (int j = 0; j < m; j++)
    {
      rez[i][j] = Determinant(GetMatr(mas, m, m, i, j), m - 1);
      if ((i + j) % 2 == 1)       
        rez[i][j] = -rez[i][j];    
      rez[i][j] = rez[i][j] / det;
    }
  }
  return Transpone(rez, m, m);
}
/*long double* cmnkB::ptlst(int km, int sn)
{
	
   long double sum;
   int i,j,k;
    
	double** c2 = new double*[sn+1];
	double* xf2 = new double[sn];
	double* xf3 = new double[sn];
	double* ca = new double[sn];

	for (int i = 0; i <= sn; i++)
	{
		c2[i] = new double[sn+1];
	}

     for (k=0;k<sn; k++)
     {
           for(j=0 ; j<sn+1; j++)
           {
                sum=0;
                for(i=0;i<km; i++) sum=sum+a[i][k]*a[i][j];
                c[k][j]=sum;
				c2[k][j]=sum;
		   }
     }
	 
	 double** c3=Mreverse(c2,sn);
	 for (k=0;k<sn; k++)
     {
		 sum=0;
		 printf("xf3[%d]=",k);
           for(j=0 ; j<sn; j++)
           {
			   sum=sum+c3[k][j]*c[j][sn];
			   //if (fabs(c3[k][j])>0.000001) printf("c3[%d][%d]=%21.18f\r\n",k,j,c3[k][j]);
			   if (fabs(c3[k][j])>0.000001) printf("c[%d]*%21.18f",j,c3[k][j]);
			   ca[j]=c[j][sn];
		   }
		   printf("\r\n");
		   //int y=0;
		   xf2[k]=sum;
	 }
	xf3[0]=ca[0]* 0.014285714285714285+ca[5]*-0.028571428571428571;
	xf3[1]=ca[1]* 0.010000000000000000;
	xf3[2]=ca[2]* 0.014285714285714285+ca[5]*-0.028571428571428571;
	xf3[3]=ca[3]* 0.020000000000000000;
	xf3[4]=ca[4]* 0.020000000000000000;
	xf3[5]=ca[0]*-0.028571428571428571+ca[2]*-0.028571428571428571+ca[5]*0.154285714285714280;
	 slvls2(sn);

	 
     
	 
	 return xf;
}*/
long double* cmnkB::ptlst(int km, int sn)
{
	
   long double sum;
   int i,j,k;

     for (k=0;k<sn; k++)
     {
           for(j=0 ; j<sn+1; j++)
           {
                sum=0;
                for(i=0;i<km; i++) sum=sum+a[i][k]*a[i][j];
                c[k][j]=sum;
		   }
     }
     slvls2(sn);
	 return xf;
}
double* cmnkB::p25()
{
	int k,i;
	double sum;
	 for (k=0;k<6; k++)
     {
         sum=0;
         for(i=0;i<25; i++) sum=sum+aa[i][k]*a25[i];
         cu[k]=sum;
     }
	xf25[0]=cu[0]* 0.014285714285714285+cu[5]*-0.028571428571428571;
	xf25[1]=cu[1]* 0.01;
	xf25[2]=cu[2]* 0.014285714285714285+cu[5]*-0.028571428571428571;
	xf25[3]=cu[3]* 0.02;
	xf25[4]=cu[4]* 0.02;
	xf25[5]=cu[0]*-0.028571428571428571+cu[2]*-0.028571428571428571+cu[5]*0.154285714285714280;
	return xf25;

}
int funccmpByte( const void * val1, const void * val2 )
{
	return *((char*) val1)>*((char*) val2);
}
int funccmpUnsignedShort( const void * val1, const void * val2 )
{
	return *((unsigned short*) val1)>*((unsigned short*) val2);
}

IEx::IEx()
{ 
	szbuf=250;
	nbufoxy=1;
	obty=0;
	aism=new int[szbuf];
	anx= new int[szbuf];
	any= new int[szbuf];
	axx=new int[szbuf];
	axy=new int[szbuf];
	aImax=new int[szbuf];
	anpix=new int[szbuf];
	ax=0;
	nstars=0;
	ax=new float[10];
	ay = new float[10];
	az = new float[10];
	zmax = new float[10];
	lnsh = new float[10];
	tetha = new float[10];
	npx =  new float[10];

	anob=0;
	omnk.initmnk(10000,6);
	nesosd=0;
	svh=0;
	svw=0;
	svt3=0;
	for (int i0=0;i0<2;i0++)
    for (int i1=0;i1<2;i1++)
	for (int i2=0;i2<2;i2++)
	for (int i3=0;i3<2;i3++)
	for (int i4=0;i4<2;i4++)
	for (int i5=0;i5<2;i5++)
	asrv[i0][i1][i2][i3][i4][i5]=0;
	
	for (int i0=0;i0<2;i0++)
    for (int i1=0;i1<2;i1++)
	for (int i2=0;i2<2;i2++)
	for (int i3=0;i3<2;i3++)
	for (int i4=0;i4<2;i4++)
	for (int i5=0;i5<2;i5++)
	for (int i6=0;i6<2;i6++)
	for (int i7=0;i7<2;i7++)
	for (int i8=0;i8<2;i8++)
	for (int i9=0;i9<2;i9++)
	asrv2[i0][i1][i2][i3][i4][i5][i6][i7][i8][i9]=0;
	char ah[5];
	char ah2[5];

	for (char i0=0;i0<5;i0++)
    for (char i1=0;i1<5;i1++)
	for (char i2=0;i2<5;i2++)
	for (char i3=0;i3<5;i3++)
	for (char i4=0;i4<5;i4++)
	{
		if ( i0!=i1&&
			 i0!=i2&&
			 i0!=i3&&
			 i0!=i4&&
			 i1!=i2&&
			 i1!=i3&&
			 i1!=i4&&
			 i2!=i3&&
			 i2!=i4&&
			 i3!=i4
			)
		{
			ah[0]=i0;
			ah[1]=i1;
			ah[2]=i2;
			ah[3]=i3;
			ah[4]=i4;
			ah2[0]=i0;
			ah2[1]=i1;
			ah2[2]=i2;
			ah2[3]=i3;
			ah2[4]=i4;
			qsort(ah2,5,1,funccmpByte);
			char srd=ah2[2];
			char ik=-1;
			for (char i5=0;i5<5;i5++)
			{
				if (ah[i5]==srd) ik=i5;
			}
			asrv2[ah[0]<ah[1]]
				 [ah[0]<ah[2]]
			     [ah[0]<ah[3]]
				 [ah[0]<ah[4]]
				 [ah[1]<ah[2]]
				 [ah[1]<ah[3]]
				 [ah[1]<ah[4]]
				 [ah[2]<ah[3]]
				 [ah[2]<ah[4]]
				 [ah[3]<ah[4]]=ik;
			int y=0;

		}
	}

	for (char i0=0;i0<4;i0++)
    for (char i1=0;i1<4;i1++)
	for (char i2=0;i2<4;i2++)
	for (char i3=0;i3<4;i3++)
	{
		if ( i0!=i1&&
			 i0!=i2&&
			 i0!=i3&&
			 i1!=i2&&
			 i1!=i3&&
			 i2!=i3
			)
		{
			ah[0]=i0;
			ah[1]=i1;
			ah[2]=i2;
			ah[3]=i3;
			
			ah2[0]=i0;
			ah2[1]=i1;
			ah2[2]=i2;
			ah2[3]=i3;
			
			qsort(ah2,3,1,funccmpByte);
			char srd1=ah2[2];
			char srd2=ah2[3];
			char ik1=-1;
			char ik2=-1;
			for (char i5=0;i5<5;i5++)
			{
				if (ah[i5]==srd1) ik1=i5;
				if (ah[i5]==srd2) ik2=i5;
			}
			asrv[ah[0]>ah[1]][ah[0]>ah[2]][ah[0]>ah[3]][ah[1]>ah[2]][ah[1]>ah[3]][ah[2]>ah[3]]=srd1;
		   asrv3[ah[0]>ah[1]][ah[0]>ah[2]][ah[0]>ah[3]][ah[1]>ah[2]][ah[1]>ah[3]][ah[2]>ah[3]]=srd2;
			

		}
	}
	//ys1=
	SetDiaFil(4);

}
void IEx::SetDiaFil(int df)
{
	t3=df;
	
	xp1=t3;
	yp1=0;
	double h=2.0*pi/5.0;
	double t=h;
	//double cost=cos(t);
	//double sint=sin(t);
	xp2=(int)(t3*cos(t)+0.5);
	yp2=(int)(t3*sin(t)+0.5);
	t=2*h;
	xp3=(int)(t3*cos(t)+0.5);
	yp3=(int)(t3*sin(t)+0.5);
	t=3*h;
	xp4=(int)(t3*cos(t)+0.5);
	yp4=(int)(t3*sin(t)+0.5);
	t=4*h;
	xp5=(int)(t3*cos(t)+0.5);
	yp5=(int)(t3*sin(t)+0.5);
	if (xp2> t3) xp2=t3;
	if (xp2<-t3) xp2=-t3;
	if (yp2> t3) yp2=t3;
	if (yp2<-t3) yp2=-t3;
	
	if (xp3> t3) xp3=t3;
	if (xp3<-t3) xp3=-t3;
	if (yp3> t3) yp3=t3;
	if (yp3<-t3) yp3=-t3;
	
	if (xp4> t3) xp4=t3;
	if (xp4<-t3) xp4=-t3;
	if (yp4> t3) yp4=t3;
	if (yp4<-t3) yp4=-t3;
	
	if (xp5> t3) xp5=t3;
	if (xp5<-t3) xp5=-t3;
	if (yp5> t3) yp5=t3;
	if (yp5<-t3) yp5=-t3;

	int y=0;
}
IEx::~IEx()
{
		delete aism;
		delete anx;
		delete any;
		delete axx;
		delete axy;
		delete aImax;
		delete anpix;
		if (anob!=0) delete anob;
		if (ax!=0)
		{ 
			delete ax;
			delete ay;
			delete az;
			delete zmax;
			delete lnsh;
			delete tetha;
			delete npx;
		}
		if (nesosd!=0)
		{
			for (int y=0;y<svh;y++)
			{
			
				delete img[y];
				delete nesosd[y];
			}
			delete img;
			delete nesosd;
			delete ap17img;
		}
		

}
void IEx::addob(int ism,int nx,int ny,int xx,int xy,int Imax,int npix)
{
	if (szbuf*nbufoxy>obty)
	{
		aism[obty]=ism;
		anx[obty]=nx;
		any[obty]=ny;
		axx[obty]=xx;
		axy[obty]=xy;
		aImax[obty]=Imax;
		anpix[obty]=npix;

		obty++;	
	} else
	{
		int* paism=aism;
		int* panx=anx;
		int* pany=any;
		int* paxx=axx;
		int* paxy=axy;
		int* paImax=aImax;
		int* panpix=anpix;
		nbufoxy++;
		aism=new int [nbufoxy*szbuf];
		anx=new int [nbufoxy*szbuf];
		any=new int [nbufoxy*szbuf];
		axx=new int [nbufoxy*szbuf];
		axy=new int [nbufoxy*szbuf];
		aImax=new int[nbufoxy*szbuf];
		anpix=new int[nbufoxy*szbuf];
		for (int i=0;i<obty;i++)
		{
			aism[i]=paism[i];
			anx[i]=panx[i];
			any[i]=pany[i];
			axx[i]=paxx[i];
			axy[i]=paxy[i];
			aImax[i]=paImax[i];
			anpix[i]=panpix[i];
		}
		delete paism;
		delete panx;
		delete pany;
		delete paxx;
		delete paxy;
		delete paImax;
		delete panpix;
		aism[obty]=ism;
		anx[obty]=nx;
		any[obty]=ny;
		axx[obty]=xx;
		axy[obty]=xy;
		aImax[obty]=Imax;
		anpix[obty]=npix;
		obty++;
	}
};

void IEx::sortob()
{
	if (anob!=0) delete anob;
	anob=new int[obty];
	for (int i=0;i<obty;i++) anob[i]=i+1;
	
	int  ii,jj;
	int x;
	
	//me:byte;
	bool hy;
	int i=obty;

     for( ii=2;ii<=i;ii++)
     {
          for( jj=i;jj>=ii;jj--)
          {

               hy=aism[jj-1-1]<aism[jj-1];
               if (hy)
			   {
                     x=aism[jj-1-1];
                     aism[jj-1-1]=aism[jj-1];
                     aism[jj-1]=x;
					 
					 x=anx[jj-1-1];
                     anx[jj-1-1]=anx[jj-1];
                     anx[jj-1]=x;
					 
					 x=any[jj-1-1];
                     any[jj-1-1]=any[jj-1];
                     any[jj-1]=x;
 					 
					 x=axx[jj-1-1];
                     axx[jj-1-1]=axx[jj-1];
                     axx[jj-1]=x;
					 
					 x=axy[jj-1-1];
                     axy[jj-1-1]=axy[jj-1];
                     axy[jj-1]=x;

 					 x=aImax[jj-1-1];
                     aImax[jj-1-1]=aImax[jj-1];
                     aImax[jj-1]=x;
					 
					 x=anob[jj-1-1];
                     anob[jj-1-1]=anob[jj-1];
                     anob[jj-1]=x;

	 				  x=anpix[jj-1-1];
                     anpix[jj-1-1]=anpix[jj-1];
                     anpix[jj-1]=x;




					 
               }
          }
	 }

};

void IEx::drawLine(int x1, int y1, int x2, int y2,int*ax,int*ay,int* n)
{
 	int ln=0;
     int deltaX = abs(x2 - x1);
     int deltaY = abs(y2 - y1);
     int signX =1;
	 if (x1 >= x2) signX = -1;
     //int signY =1; y1 < y2 ? 1 : -1;
	 int signY =1;
	 if (y1 >= y2) signY = -1;
    //
    int error = deltaX - deltaY;
    //
    //setPixel(x2, y2);
	int sx2=x2;
	int sy2=y2;
	
    while(x1 != x2 || y1 != y2) 
   {
        //setPixel(x1, y1);
		ax[ln]=x1;
		ay[ln]=y1;
		ln++;
        const int error2 = error * 2;
        //
        if(error2 > -deltaY) 
        {
            error -= deltaY;
            x1 += signX;
        }
        if(error2 < deltaX) 
        {
            error += deltaX;
            y1 += signY;
        }
    }
	ax[ln]=sx2;
	ay[ln]=sy2;
	ln++;
	*n=ln;

}
int funccmp( const void * val1, const void * val2 )
{
	return *((double*) val1)>*((double*) val2);
}
void IEx::GetCen2(unsigned short** ap17img,unsigned short** nesosdQ,int tt100, int w, int h,int t765,int t65)
{
	tt100=obty;
	if (ax!=0)
	{ 
			delete ax;
			delete ay;
			delete az;
			delete zmax;
			delete lnsh;
			delete tetha;
			delete npx;
	}
	ax=new float[tt100];
	ay=new float[tt100];
	az=new float[tt100];
	zmax=new float[tt100];
	lnsh=new float[tt100];
	tetha=new float[tt100];
	npx=new float [tt100];


	nstars=0;
	
	
	//tt100=100;
	int dxa[8]={-1, 0, 1,-1,1,-1,0,1};
	int dya[8]={-1,-1,-1, 0,0, 1,1,1};


	int nnas=0;
	int nno=0;
	int i=0;
	int ngo=0;
	//for (int i=tt100-200;i<tt100;i++)
	//for (int i=0;i<tt100;i++)
	//tt100=1;
	while (nno<t65&&
		i<tt100)
	{
		int nob=anob[i];
		//aImax[i]=0;
		int bx=anx[i];
		int ex=axx[i]+1;
		int by=any[i];
		int ey=axy[i]+1;
		int nu=0;
		int max=0;
		int mx=-1;
		int my=-1;
		bool shtrih=false;
		if (anpix[i]>16)
		{
	
			int dxb=ex-bx;
			int dyb=ey-by;
			int sb=dxb*dyb;
		
			bool** nesosdW=new bool*[dyb];
			for (int y=0;y<dyb;y++)
			{
				nesosdW[y]=new bool[dxb];
				for (int x=0;x<dxb; x++)
				{
					nesosdW[y][x]=false;
				}
			}
			int nu=0;
			//for (int y=by;y<ey;y++)
			{
			
				int y=by;
				for (int x=bx;x<ex;x++)
				{
					if (nesosdQ[y][x]==nob)
					{
						nesosdW[y-by][x-bx]=true;
						nu++;
					}
				}
				y=ey-1;
				for (int x=bx;x<ex;x++)
				{
					if (nesosdQ[y][x]==nob)
					{
						nesosdW[y-by][x-bx]=true;
						nu++;
					}
				}
				
			}
			{
			
				int x=bx;
				for (int y=by;y<ey;y++)
				{
					if (nesosdQ[y][x]==nob)
					{
						nesosdW[y-by][x-bx]=true;
						nu++;
					}
				}
				x=ex-1;
				for (int y=by;y<ey;y++)
				{
					if (nesosdQ[y][x]==nob)
					{
						nesosdW[y-by][x-bx]=true;
						nu++;
					}
				}
				
			}
			for (int y=by+1;y<ey-1;y++)
			{
				for (int x=bx+1;x<ex-1;x++)
				{
					
					if (nesosdQ[y][x]==nob)
					{
						bool bgr=false;
						for (int i=0;i<8;i++) 
						{
							bgr=bgr||nesosdQ[y+dya[i]][x+dxa[i]]!=nob;
						}
						if (bgr)
						{
							nesosdW[y-by][x-bx]=true;
							nu++;
						}
					}
				}
			}
			
			double sx=0;
			double sy=0;
			int nuu=0;
			for (int y=by;y<ey;y++)
			{
				for (int x=bx;x<ex;x++)
				{
					if (nesosdW[y-by][x-bx])
					{
					   sx=sx+x;
					   sy=sy+y;
					   nuu++;
					   
					}
				}
			}
			if (nuu>0)
			{
				sx=sx/nuu;
				sy=sy/nuu;
			}
			double maxr=0;
			double minr=1e39;
			int mxx=0;
			int mxy=0;
			for (int y=by;y<ey;y++)
			{
				for (int x=bx;x<ex;x++)
				{
					if (nesosdW[y-by][x-bx])
					{
						double dx=sx-x;
						double dy=sy-y;
						double r=dx*dx+dy*dy;
						if (r>maxr)
						{
							maxr=r;
							mxx=x;
							mxy=y;
						}
						if (r<minr) minr=r;

						
					   
					}
				}
			}
			minr=sqrt(minr);
			maxr=sqrt(maxr);
			if (minr<1) minr=1;
			double otn=maxr/minr;
			
			/*{
			for (int y=by;y<ey;y++)
			{
				for (int x=bx;x<ex;x++)
				{
					if (nesosdW[y-by][x-bx]) ap17img[y][x]=0;
				}
			}
			}*/
			if (otn>4 && maxr>4)
			{
				shtrih=true;
				int* axL=new int[((int)maxr)*10];
				int* ayL=new int[((int)maxr)*10];

				int isx=(int)(sx+0.5);
				int isy=(int)(sy+0.5);
				maxr=0;
				int mx1s=-1;
				int  mx2s=-1;
				int  my1s=-1;
				int  my2s=-1;

				for (int y=by;y<ey;y++)
				{
					for (int x=bx;x<ex;x++)
					{
						if (nesosdW[y-by][x-bx])
						{
							int dx=isx-x;
							int dy=isy-y;
							int x1=isx-dx*2;
							int y1=isy-dy*2;
							int x2=isx+dx*2;
							int y2=isy+dy*2;
							int n=0;
							drawLine(x1,y1,x2,y2,axL,ayL,&n);
							int x1s=-1;
							int y1s=-1;
							int m=0;
							do
							{
								if (
									bx<=axL[m]&&axL[m]<ex &&
									by<=ayL[m]&&ayL[m]<ey
									)
								{
									if (nesosdW[y-by][x-bx])
									{
										x1s=axL[m];
										y1s=ayL[m];
									}

								}
								m++;
							} while (m<n&&x1s<0);
							
							int x2s=-1;
							int y2s=-1;
							m=n-1;
							do
							{
								if (
									bx<=axL[m]&&axL[m]<ex &&
									by<=ayL[m]&&ayL[m]<ey
									)
								{
									if (nesosdW[y-by][x-bx])
									{
										x2s=axL[m];
										y2s=ayL[m];
									}

								}
								m--;
							} while (m>=0&&x2s<0);
							if (x1s>=0&&x2s>=0)
							{
								dx=x2s-x1s;
								dy=y2s-y1s;
								double r=dx*dx+dy*dy;
								if (r>maxr)
								{
									maxr=r;
									mx1s=x1s;
									mx2s=x2s;
									my1s=y1s;
									my2s=y2s;
								}

							}
							int ys=0;
						}
					}
				}
				int n=0;
				//drawLine(mx1s,my1s,mx2s,my2s,ax,ay,&n);
				/*double x0s=(mx1s+mx2s)/2.0;
				double y0s=(my1s+my2s)/2.0;
				FILE* file = fopen("otl.txt", "wb");
				for (int i=0;i<n;i++)
				{
					//ap17img[ay[i]][ax[i]]=0;
					double dx=(ax[i]-mx1s);
					double dy=(ay[i]-my1s);
					double r=sqrt(dx*dx+dy*dy);
					//CString pst;
					//pst.Format("%10.5f %5d\r\n",r,ap17img[ay[i]][ax[i]]);
					fprintf(file,"%10.5f %5d %5d %5d\r\n",r,ap17img[ay[i]][ax[i]],ax[i],ay[i]);
				}
				fclose(file);*/
				double* ar=new double [nuu];
				double dx=(mx2s-mx1s);
				double dy=(my2s-my1s);
				double r=sqrt(dx*dx+dy*dy);
				double u1=mx2s*my1s-my2s*mx1s;
				int nt=0;
				for (int y=by;y<ey;y++)
				{
					for (int x=bx;x<ex;x++)
					{
						if (nesosdW[y-by][x-bx])
						{
							
							ar[nt]=fabs(x*dy-y*dx+u1)/r;
							//ap17img[y][x]=0;
							nt++;

						}
					}
				}
				qsort(ar,nt,8,funccmp);
				double srd=ar[nt/2]*1.8;
				if (srd<2.85) srd=2.85;
				//srd=15;
				//srd=srd+1;
				
				int dr=((int)srd)+1;
				int bx2=bx-dr-1;
				int ex2=ex+dr+2;
				int by2=by-dr-1;
				int ey2=ey+dr+2;
				
				if (by2<0) by2=0;
				if (by2>h) by2=h;

				if (ey2<0) ey2=0;
				if (ey2>h) ey2=h;
				if (bx2<0) bx2=0;
				if (bx2>w) bx2=w;
				if (ex2<0) ex2=0;
				if (ex2>w) ex2=w;

				int x1=mx1s;
				int x2=mx2s;
				int y1=my1s;
				int y2=my2s;
				int x1a=x1;
				int x2a=x2;
				int y1a=y1;
				int y2a=y2;
				if (x1>x2)
				{
					x1a=x2;
					x2a=x1;
				}
				if (y1>y2)
				{
					y1a=y2;
					y2a=y1;
				}
				double a=y1-y2;
				double b=x2-x1;
				double c=x1*y2-x2*y1;
				
				double bb=b*b;
				double ac=a*c;
				double bc=b*c;
				double ab=a*b;
				double aa=a*a;
				double zab=aa+bb;
				srd=srd*srd;
				int km=0;
				for (int y=by2;y<ey2;y++)
				{
					int y2=y-by2;
					for (int x=bx2;x<ex2;x++)
					{
						bool bur=false;
						double xt= (bb*x-ab*y-ac)/zab;
						double yt=(-ab*x+aa*y-bc)/zab;
						if ( x1a<=xt&&xt<=x2a&&y1a<=yt&&yt<=y2a)
						{
							double r1=fabs(x*dy-y*dx+u1)/r;
							r1=r1*r1;
							if (r1<srd)
							{
								//ap17img[y][x]=0;
								bur=true;
							}
							

						} else
						{
							double dx2=x-mx1s;
							double dy2=y-my1s;
							double dx3=x-mx2s;
							double dy3=y-my2s;
							double r2=dx2*dx2+dy2*dy2;
							double r3=dx3*dx3+dy3*dy3;
							if (r3<r2) r2=r3;
							if (r2<srd)
							{
								//ap17img[y][x]=0;
								bur=true;
							}
							
						}
						if (km<10000&& bur&&ap17img[y][x]<t765)
						{
							int x2=x-bx2;
							omnk.a[km][0]=x2*x2;
							omnk.a[km][1]=x2*y2;
							omnk.a[km][2]=y2*y2;
							omnk.a[km][3]=x2;
							omnk.a[km][4]=y2;
							omnk.a[km][5]=1;
							omnk.a[km][6]=ap17img[y][x];
							km++;
						}
						
						
					}
					
				}
				if (km<10000&&km>6)
				{
					long double* b=omnk.ptlst(km,6);
					if (b[1]!=0)
					{
						double c1=-4*b[0]*b[2]/b[1];
						if (c1+b[1]!=0)
						{
							double c2=-(b[4]/b[1])*b[0]*2+b[3];
							double yc=-c2/(c1+b[1]);
							double xc=-(2*b[2]*yc+b[4])/b[1]+bx2;
							yc=yc+by2;
							//ap17img[(int)yc][(int)xc]=0;
							double dxw=x1-xc;
							double dyw=y1-yc;
							double dfte=sqrt(dyw*dyw+dxw*dxw);
							double tetaa=atan(dxw/dyw);
							if (dyw<0) tetaa=pi+tetaa;
							if (tetaa<0) tetaa=tetaa+2*pi;
							tetaa=tetaa*180/pi;
							double dte=180/(dfte*pi);
							double dte0=dte;
							int itr=0;
							double sm2=1;
							double bu[6];
							for (int q=0;q<6;q++) bu[q]=b[q];
							//tetaa=0;
							bool vyh=false;

							do
							{
								km=0;
								for (int du=-2;du<3;du++)
								{
									double te=tetaa+du*dte;
									double x=dfte*sin(te*pi/180)+xc-bx2;
									double y=dfte*cos(te*pi/180)+yc-by2;
									double z=bu[0]*x*x+bu[1]*x*y+bu[2]*y*y+x*bu[3]+y*bu[4]+bu[5];
									omnk.a[km][0]=te*te;
									omnk.a[km][1]=te;
									omnk.a[km][2]=1;
									omnk.a[km][3]=z;
									km++;
								}
								long double * bq=omnk.ptlst(km,3);
								long double * sgm=omnk.sgmnk(km,3);
								double tn=0;
								if (bq[0]!=0)
								{
								   tn=-bq[1]/(2*bq[0]);
								} else
								{
									vyh=true;
								}
								double sm1=fabs(tn-tetaa)/dte;
								if (sm1<2) dte=dte/2;
								sm2=fabs(tn-tetaa)/dte0;
								tetaa=tn;
								itr++;
							} while (sm2>0.001 && itr<10000&&!vyh);
							if (itr<10000 &&!vyh)
							{
								double xc2=dfte*sin(tetaa*pi/180)+xc;
								double yc2=dfte*cos(tetaa*pi/180)+yc;
								double a=yc-yc2;
								double b=xc2-xc;
								double c=xc*yc2-xc2*yc;
				
								double bb=b*b;
								double ac=a*c;
								double bc=b*c;
								double ab=a*b;
								double aa=a*a;
								double zab=aa+bb;
								if (zab!=0)
								{
									double xt1= (bb*x1-ab*y1-ac)/zab;
									double yt1=(-ab*x1+aa*y1-bc)/zab;
									
									double xt2= (bb*x2-ab*y2-ac)/zab;
									double yt2=(-ab*x2+aa*y2-bc)/zab;
									
									double xc2=(xt1+xt2)*0.5;
									double yc2=(yt1+yt2)*0.5;
									
									double dx=xt1-xt2;
									double dy=yt1-yt2;
									double r=sqrt(dx*dx+dy*dy);
									ax[nstars]=xc2;
									ay[nstars]=yc2;
									az[nstars]=aism[i];
									zmax[nstars]=aImax[i];
									lnsh[nstars]=r;
									tetha[nstars]=tetaa;
									npx[nstars]=anpix[i];
									nstars++;
								}
							}

							
						}
					}
					/*for (int y=by2;y<ey2;y++)
					{
						int y2=y-by2;
						for (int x=bx2;x<ex2;x++)
						{
							bool bur=false;
							double xt= (bb*x-ab*y-ac)/zab;
							double yt=(-ab*x+aa*y-bc)/zab;
							if ( x1a<=xt&&xt<=x2a&&y1a<=yt&&yt<=y2a)
							{
								double r1=fabs(x*dy-y*dx+u1)/r;
								r1=r1*r1;
								if (r1<srd)
								{
									//ap17img[y][x]=0;
									bur=true;
								}
							

							} else
							{
								double dx2=x-mx1s;
								double dy2=y-my1s;
								double dx3=x-mx2s;
								double dy3=y-my2s;
								double r2=dx2*dx2+dy2*dy2;
								double r3=dx3*dx3+dy3*dy3;
								if (r3<r2) r2=r3;
								if (r2<srd)
								{
									//ap17img[y][x]=0;
									bur=true;
								}
							
							}
							if (km<10000&& bur&&ap17img[y][x]<t765)
							{

								int x2=x-bx2;
								long double* xf=b;
								double z=x2*x2*xf[0]+x2*y2*xf[1]+y2*y2*xf[2]+x2*xf[3]+y2*xf[4]+xf[5];
								ap17img[y][x]=z;
								
							}
						
						
						}
					
					}*/
				}
				
				//ap17img[my2s][mx2s]=765;

				

				delete ar;
				delete axL;
				delete ayL;

			}
			


			
			for (int y=0;y<ey-by;y++)
			{
				delete nesosdW[y];
				
			}
			delete nesosdW;
			//return;
		}
		if (!shtrih)
		{
		
		

		for (int y=by;y<ey;y++)
		{
			for (int x=bx;x<ex;x++)
			{
				if (nesosdQ[y][x]==nob)
				{
					if (ap17img[y][x]>max)
					{
						max=ap17img[y][x];
						mx=x;
						my=y;
						nu++;
					}
				}
				
			}
		}
		//if (max==t765) nnas++;
		
		
		if (max<t765)
		{
			int km=0;
			for (int y3=-2;y3<3;y3++)
			{
				for (int x3=-2;x3<3;x3++)
				{
					/*omnk.a[km][0]=x3*x3;
					omnk.a[km][1]=x3*y3;
					omnk.a[km][2]=y3*y3;
					omnk.a[km][3]=x3;
					omnk.a[km][4]=y3;
					omnk.a[km][5]=1;
					omnk.a[km][6]=ap17img[my+y3][mx+x3];*/
					omnk.a25[km]=ap17img[my+y3][mx+x3];
					km++;
				}
			
			}
			//long double* xfx=omnk.ptlst(km,6);
			//long double* sgm=omnk.sgmnk(km,6);
			double* b=omnk.p25();
			double dy=10000;
			double dx=10000;
				if (b[1]!=0)
				{
					
					double c1=-4*b[0]*b[2]/b[1];
					if (c1+b[1]!=0)
					{
						double c2=-(b[4]/b[1])*b[0]*2+b[3];
						dy=-c2/(c1+b[1]);
						dx=-(2*b[2]*dy+b[4])/b[1];
						int y=0;
					}
				}
			
			//double px=2*b[0]*dx+b[1]*dy+b[3];
			//double py=b[1]*dx+2*b[2]*dy+b[4];
			double r2=dx*dx+dy*dy;
			if (r2<4
				//||true
				)
			{
									ax[nstars]=mx+dx;
									ay[nstars]=my+dy;
									az[nstars]=aism[i];
									zmax[nstars]=aImax[i];
									lnsh[nstars]=0;
									tetha[nstars]=0;
									npx[nstars]=anpix[i];
									nstars++;

				ngo++;
				/*for (int y=by;y<ey;y++)
				{
					for (int x=bx;x<ex;x++)
					{
				 		if (nesosdQ[y][x]==nob)
						{
							ap17img[y][x]=0;
						}
					}
				}*/
				
			

			
			} else
			{
				
				nno++;
			}
		} else
		{
			int bx=anx[i];
			int ex=axx[i]+1;
			int by=any[i];
			int ey=axy[i]+1;
			
			int dx=ex-bx+6;
			int dy=ey-by+6;
			
			int bxb=bx-3;
			int exb=ex+3;
			int byb=by-3;
			int eyb=ey+3;
			bool** bobl=new bool*[dy];
			for (int y=0;y<dy;y++)
			{
				bobl[y]=new bool[dx];
				for (int x=0;x<dx;x++)
				{
					bobl[y][x]=false;
				}
			}
			for (int y=by;y<ey;y++)
			{
				for (int x=bx;x<ex;x++)
				{
					if (ap17img[y][x]>=t765)
					{
						int bgra=false;
						for (int i=0;i<8;i++)
						{
							bgra=bgra||ap17img[y+dya[i]][x+dxa[i]]<t765;
						}
						if (bgra)
						{
							for (int y2=y-2;y2<y+3;y2++)
							{
								for (int x2=x-2;x2<x+3;x2++)
								{
									bobl[y2-byb][x2-bxb]=bobl[y2-byb][x2-bxb]||ap17img[y2][x2]<t765;
									int yd=0;

								}
							}

						}
					}

				}
			}
			int km=0;
			for (int y=byb;y<eyb;y++)
			{
				int y2=y-byb;
				for (int x=bxb;x<exb;x++)
				{
					if (bobl[y-byb][x-bxb]&&km<10000)
					{
						int x2=x-bxb;
						omnk.a[km][0]=x2*x2;
						omnk.a[km][1]=x2*y2;
						omnk.a[km][2]=y2*y2;
						omnk.a[km][3]=x2;
						omnk.a[km][4]=y2;
						omnk.a[km][5]=1;
						omnk.a[km][6]=ap17img[y][x];
						km++;
					}
					
					
				}
			
			}
			if (km<10000&&km>6)
			{
				long double* b=omnk.ptlst(km,6);
				//long double* sgm=omnk.sgmnk(km,6);
				if (b[1]!=0)
				{
					double c1=-4*b[0]*b[2]/b[1];
					if (c1+b[1]!=0)
					{
						double c2=-(b[4]/b[1])*b[0]*2+b[3];
						double yc=-c2/(c1+b[1]);
						double xc=-(2*b[2]*yc+b[4])/b[1]+bxb;
						yc=yc+byb;
									ax[nstars]=xc;
									ay[nstars]=yc;
									az[nstars]=aism[i];
									zmax[nstars]=aImax[i];
									lnsh[nstars]=0;
									tetha[nstars]=0;
									npx[nstars]=anpix[i];
									nstars++;

						int o=0;
					}
				}
			}
			/*		
			for (int y=byb;y<eyb;y++)
			{
				for (int x=bxb;x<exb;x++)
				{
					if (bobl[y-byb][x-bxb])
					{
						ap17img[y][x]=0;

					}
				}
			}*/

			for (int y=0;y<dy;y++)
			{
				delete bobl[y];
			}
			delete bobl;
			


			
		}
		}
				/*for (int y=by;y<ey;y++)
				{
					for (int x=bx;x<ex;x++)
					{
				 		if (nesosdQ[y][x]==nob)
						{
							
							ap17img[y][x]=0;
						}
					}
				}*/
		i++;

	}
	int y=0;
}
void OneSour(int* ism,int* nx,int* ny,int* xx,int* xy, int x,int y,long** img3,unsigned short** nesosd,int bx,int by,int ex,int ey,int z,int* glubrek,int* Imax,unsigned short**ap17img,int* mxi,int*myi,int urnas,int nob,int* npix)
{
	//if (glubrek>6425)
	//{
	//	int y=0;
	//}
	//if (glubrek>mxgrek) mxgrek=glubrek;

	(*glubrek)++;
	if (!nesosd[y][x])
	{ 
		int dx[8]={-1, 0, 1,-1,1,-1,0,1};
		int dy[8]={-1,-1,-1, 0,0, 1,1,1};
		nesosd[y][x]=nob+1;
		(*npix)++;
		*ism=*ism+img3[y][x];
		if (x<*nx) *nx=x;
		if (x>*xx) *xx=x;
		if (y<*ny) *ny=y;
		if (y>*xy) *xy=y;
		bool uge=false;
		if (ap17img[y][x]>*Imax)
		{
			*Imax=ap17img[y][x];
			*mxi=x;
			*myi=y;
		}
		do
		{
			uge=false;
			int sx=x;
			int sy=y;
			for (int i=0;i<8;i++)
			{
				int y2=sy+dy[i];
				int x2=sx+dx[i];
				if (bx<=x2&&x2<=ex&&by<=y2&&y2<=ey)
				{
					if ((img3[y2][x2]>z ||ap17img[y2][x2]>=urnas)&& !nesosd[y2][x2])
					{
						if (uge)
						{
							if ((*glubrek)<1000)
							{
								OneSour(ism,nx,ny,xx,xy,x2,y2,img3,nesosd,bx,by,ex,ey,z,glubrek,Imax,ap17img,mxi,myi,urnas,nob,npix);
							}
						} else
						{
							x=x2;
							y=y2;
							nesosd[y2][x2]=nob+1;
							(*npix)++;
							*ism=*ism+img3[y2][x2];
							if (img3[y2][x2]>*Imax)
							{
								*Imax=img3[y2][x2];
								*mxi=x2;
								*myi=y2;
							}


							if (x2<*nx) *nx=x2;
							if (x2>*xx) *xx=x2;
							if (y2<*ny) *ny=y2;
							if (y2>*xy) *xy=y2;
							uge=true;
						}
			
					}
				}
			}
		} while (uge);
	
	}
	(*glubrek)--;
}
void qsortRecursive(int *mas, int size) 
{
    
    int i = 0;
    int j = size - 1;

    
    int mid = mas[size / 2];

    
    do {
        while(mas[i] < mid) {
            i++;
        }
    
        while(mas[j] > mid) {
            j--;
        }

        
        if (i <= j) {
            int tmp = mas[i];
            mas[i] = mas[j];
            mas[j] = tmp;

            i++;
            j--;
        }
    } while (i <= j);


    
    if(j > 0) {
    
        qsortRecursive(mas, j + 1);
    }
    if (i < size) {
    
        qsortRecursive(&mas[i], size - i);
    }
}

void FusPix(long** img3,unsigned short** nesosd,int bx,int by,int ex,int ey,int z,IEx* ccd1,unsigned short**ap17img,int urnas)
{
	int dxr[]={-3,-3,-2,-2,-1,-1, 0,0, 1,1, 2,2,3};
	int dyr[]={ 1, 0,-1, 2,-1, 2,-2,3,-1,2,-1,2,0};
	int zz[13];

	for(int y=by;y<=ey;y++)
	{
		for (int x=bx;x<=ex;x++)
		{
			if (y==768&&x==1199)
			{
				int yf=0;
			}

			if ( (img3[y][x]>z||ap17img[y][x]>=urnas)&&!nesosd[y][x])
			{
				int nx,ny,xx,xy,ism,Imax;
				nx=1000000;
				ny=1000000;
				xx=-1000000;
				xy=-1000000;
				Imax=-1000000;
				ism=0;
				int glubrek=0;
				int mxi,myi;
				int npix=0;
				OneSour(&ism,&nx,&ny,&xx,&xy,x,y,img3,nesosd,bx,by,ex,ey,z,&glubrek,&Imax,ap17img,&mxi,&myi,urnas,ccd1->obty,&npix);
				for (int i=0;i<13;i++)
				{
					zz[i]=ap17img[myi+dyr[i]][mxi+dxr[i]];
				}
				qsortRecursive(zz,13);
				Imax=Imax-zz[6];
				
				ccd1->addob(ism,nx,ny,xx,xy,Imax,npix);
				

			}
		}
	}
	//mxgrek++;
}
long** IEx::GetSubBG()
{
	return img;
}




void IEx::GetSources(unsigned short* ccdimgC, int w, int h,double dS,int nLimSh ,int urnas,
	int* nstr,
	float** axS,
	float**	ayS,
	float** azS,
	float** zmaxS,
	float** lnshS,
	float** tethaS,
	float** npxS)


{
	
	
	int t4g=4;
	int tk=3;
	//int urnas=255;
	int t400=(int)w*h*dS;//
	

	//ap17img[y+dy][x+dx];
	int dxr[]={-3,-3,-2,-2,-1,-1, 0,0, 1,1, 2,2,3};
	int dyr[]={ 1, 0,-1, 2,-1, 2,-2,3,-1,2,-1,2,0};
	//int zz[13];
	obty=0;
	nstars=0;
	if (nesosd==0 || h!=svh|| w!=svw|| t3!=svt3)
	{
		if (nesosd!=0)
		{
			for (int y=0;y<svh;y++)
			{
			
				delete img[y];
				delete nesosd[y];
			}
			delete img;
			delete nesosd;
			delete ap17img;
		}

	
		nesosd=new unsigned short* [h];

		for(int i=0; i<h;i++)
		{
			nesosd[i]=new unsigned short[w];
			for (int j=0;j<w;j++)
			{
				nesosd[i][j]=0;
			}

		}
	
	
		ap17img= new unsigned short*[h];
		img= new long*[h];
		for (int y=0;y<h;y++)
		{
			img[y]=new long[w];
		}
	
		for (int y=0;y<t3;y++)
		{
			for (int x=0;x<w;x++)
			{
				img[y][x]=0;
			}
		}
		for (int y=h-t3;y<h;y++)
		{
			for (int x=0;x<w;x++)
			{
				img[y][x]=0;
			}
		}

		for (int x=0;x<t3;x++)
		{
		
			for (int y=0;y<h;y++)
			{
				img[y][x]=0;
			}
		}
	
		for (int x=w-t3;x<w;x++)
		{
			for (int y=0;y<h;y++)	
			{
				img[y][x]=0;
			}
		}
		svh=h;
		svw=w;
		svt3=t3;
	}
	for (int y=0;y<h;y++)
	{
		ap17img[y]=ccdimgC+y*w;
	}

	for(int i=0; i<h;i++)
	{
		
			for (int j=0;j<w;j++)
			{
				nesosd[i][j]=0;
			}

	}
	
	//int z1,z2,z3,z4,max,min;
	//unsigned short  a,b,c,d;
	//unsigned short  swap;
	//long  a,b,c,d;
	//int  swap;
	unsigned short ah[5];
	long a;
	//int u2=0;
	char axp[6];
	short ayt[6];
	axp[1]=xp1;
	axp[2]=xp2;
	axp[3]=xp3;
	axp[4]=xp4;
	axp[5]=xp5;
	for (int y=t3;y<h-t3;y=y+2)
	{
		
		short yt1=y+yp1;
		short yt2=y+yp2;
		short yt3=y+yp3;
		short yt4=y+yp4;
		short yt5=y+yp5;

		for (int x=t3;x<w-t3;x=x+2)
		{
				
				ah[0]=ap17img[yt1][x+axp[1]];
				ah[1]=ap17img[yt2][x+axp[2]];
				ah[2]=ap17img[yt3][x+axp[3]];
				ah[3]=ap17img[yt4][x+axp[4]];
				ah[4]=ap17img[yt5][x+axp[5]];
				
				a=ah[asrv2[ah[0]<ah[1]]
				          [ah[0]<ah[2]]
						  [ah[0]<ah[3]]
						  [ah[0]<ah[4]]
						  [ah[1]<ah[2]]
						  [ah[1]<ah[3]]
						  [ah[1]<ah[4]]
						  [ah[2]<ah[3]]
						  [ah[2]<ah[4]]
						  [ah[3]<ah[4]]
				     ];
				
				img[y][x]=ap17img[y][x]-a;
				img[y][x+1]  =ap17img[y  ][x+1]-a;
				img[y+1][x]  =ap17img[y+1][x ]-a;
				img[y+1][x+1]=ap17img[y+1][x+1]-a;
		}
	}
	int max=img[0][0];
	int min=img[0][0];
	for (int y=0;y<h;y=y+t4g)
	{
		for (int x=0;x<w;x=x+t4g)
		{
			if (img[y][x]>max) max=img[y][x];
			if (img[y][x]<min) min=img[y][x];
		}
	}


	int dmm=max-min+1;
	unsigned short* gist=new unsigned short [dmm];
	for (int ie=0;ie<dmm;ie++)
	{
		gist[ie]=0;
	}
	for (int y=0;y<h;y=y+t4g)
	{
		for (int x=0;x<w;x=x+t4g)
		{
			int z=img[y][x]-min;
			gist[z]=gist[z]+1;
		}
	}
	int z=max;
	int ns=0;
	while (ns<t400 && z>min)
	{
		
		int z2=z-min;
		ns=ns+gist[z2];
		z--;
	}
	if (ns>t400*2)
	{
		z++; 
	}
	delete gist;
	//IEx ccd1;
	
	FusPix(img,nesosd,t3,t3,w-t3,h-t3,z,this,ap17img,urnas);
	sortob();
	//*numk=anpix[0];
	GetCen2(ap17img,nesosd,0,w, h,urnas, nLimSh);
	*nstr = nstars;
	*axS = ax;
	*ayS = ay;
	*azS = az;
	*zmaxS = zmax;
	*lnshS = lnsh;
	*tethaS = tetha;
	*npxS = npx;
	
	

	/*int z=300;
	int nz=0;
	for (int y=3;y<h-3;y++)
	{
		for (int x=3;x<w-3;x++)
		{
			if (img[y][x]>z) nz++;
		}
	}*/
	//if (true)
	//if (false)
	/*if (ccd1.anpix[0]>100)
	{
		FILE* file = fopen(pFileName, "wb");

		if (file)
		{
      // RGB image
			int imageSizeInBytes = 3 * w * h;
			int headersSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
			int fileSize = headersSize + imageSizeInBytes;
			
			uint8_t * pData = new uint8_t[headersSize];
			
			if (pData != NULL)
			{
				BITMAPFILEHEADER& bfHeader = *((BITMAPFILEHEADER *)(pData));

				bfHeader.bfType = 'MB';
				bfHeader.bfSize = fileSize;
				bfHeader.bfOffBits = headersSize;
				bfHeader.bfReserved1 = bfHeader.bfReserved2 = 0;

				BITMAPINFOHEADER& bmiHeader = *((BITMAPINFOHEADER *)(pData + headersSize - sizeof(BITMAPINFOHEADER)));

				bmiHeader.biBitCount = 3 * 8;
				bmiHeader.biWidth    = w;
				bmiHeader.biHeight   = h;
				bmiHeader.biPlanes   = 1;
				bmiHeader.biSize     = sizeof(bmiHeader);
				bmiHeader.biCompression = BI_RGB;
				bmiHeader.biClrImportant = bmiHeader.biClrUsed = 
          bmiHeader.biSizeImage = bmiHeader.biXPelsPerMeter = 
          bmiHeader.biYPelsPerMeter = 0;
				
        fwrite(pData, headersSize, 1, file);
		int numberOfBytesToWrite = 3 * w;

				unsigned char *pBits = new unsigned char [numberOfBytesToWrite];
				int nSpan      = - frame->linesize[0];
				

				

				for (size_t i = 0; i < h; i++)	
				{
					for (int j=0;j<w;j++)
					{
						//int z2=img[i][j]/3+100;
						int z2=ap17img[i][j]/3;
						if (z2<0) z2=0;
						if (z2>255) z2=255;
						//if (img[i][j]<z) z2=0;
						pBits[j*3]=z2;
						pBits[j*3+1]=z2;
						pBits[j*3+2]=z2;
					}
					fwrite(pBits, numberOfBytesToWrite, 1, file);
				}
				delete pBits;
				
				
				delete [] pData;				
			}

			fclose(file);
		}
	}*/

	
}



extern "C"
{
    void szm(unsigned short* ccdimgC, int w, int h, float* x, float* y, int* nstr)
    {
        IEx OEx;
        float* ax;
        float* ay;
        float* az;
        float* zmax;
        float* lnsh;
        float* tetha;
        float* npx;
        OEx.GetSources(ccdimgC,//input buffer (ccd-frame)
            w,h,
            0.0002,//fraction source-occuped square
            1000,//accepted number of false detections
            20000,//saturation threshold value
            nstr,// number of detected sources
            &ax,//X-pisition output array
            &ay,//Y-pisition output array
            &az,// flux array
            &zmax,//maximum signal value
            &lnsh,//segment lenght (if sourse classified as a line segment)
            &tetha,//segment positional angle
            &npx// number of pixels
            );
        int ns = nstr[0];
        for(int k=0;k<ns;k++){x[k]=ax[k];y[k]=ay[k];}
    }
    
    void szm2(unsigned short* ccdimgC, int w, int h, float* x, float* y, float* l, float* angs, int* nstr)
    {
        IEx OEx;
        //int nstr;
        // input arrays
        float* ax;
        float* ay;
        float* az;
        float* zmax;
        float* lnsh;
        float* tetha;
        float* npx;
        OEx.GetSources(ccdimgC,//input buffer (ccd-frame)
            w,h,
            0.0003,//fraction source-occuped square
            5000,//accepted number of false detections
            20000,//saturation threshold value
            nstr,// number of detected sources
            &ax,//X-pisition output array
            &ay,//Y-pisition output array
            &az,// flux array
            &zmax,//maximum signal value
            &lnsh,//segment lenght (if sourse classified as a line segment)
            &tetha,//segment positional angle
            &npx// number of pixels
            );
        //copy to output arrays
        int ns = nstr[0];
        for(int k=0;k<ns;k++){x[k]=ax[k];y[k]=ay[k];l[k]=lnsh[k];angs[k]=tetha[k];}
    }
    
}
