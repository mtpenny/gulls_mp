#include<cmath>
#include<vector>

#include "constants.h"
#include "coords.h"
#include "integerPowers.h"

using namespace std;

void coords::ad2lb(double a, double d, double* l, double* b)
{
  vector<double> lb;
  ad2lb(a,d,&lb);
  *l = lb[0]; *b = lb[1];
}

void coords::ad2lb(double a, double d, vector<double>* lb)
{
  double cl, sl, cb, sb, cd, sd, ca, sa;
  double xgal, ygal;
  lb->resize(6);

  ca=cos(a); sa=sin(a);
  cd=cos(d); sd=sin(d);

  xgal = T[0][0]*cd*ca + T[0][1]*cd*sa + T[0][2]*sd;
  ygal = T[1][0]*cd*ca + T[1][1]*cd*sa + T[1][2]*sd;
  sb = T[2][0]*cd*ca + T[2][1]*cd*sa + T[2][2]*sd;
  cb = hypot(xgal, ygal);
  (*lb)[1] = atan2(sb,cb);
  //(*lb)[1] = asin(sb);
  //cb = cos((*lb)[1]);

  cl = xgal/cb;
  sl = ygal/cb;
  (*lb)[0] = atan2(sl,cl);
  if((*lb)[0]<0){(*lb)[0]+=twoPi;}
  if((*lb)[0]>twoPi){(*lb)[0]-=twoPi;}

  (*lb)[2]=cl; (*lb)[3]=sl;   (*lb)[4]=cb; (*lb)[5]=sb;
}

void coords::lb2ad(double l, double b, double* a, double* d)
{
  vector<double> ad;
  lb2ad(l,b,&ad);
  *a = ad[0]; *d = ad[1];
}

void coords::lb2ad(double l, double b, vector<double>* ad)
{
  double cl, sl, cb, sb, cd, sd, ca, sa;
  double xeq, yeq;
  ad->resize(6);

  cl=cos(l); sl=sin(l);
  cb=cos(b); sb=sin(b);

  xeq = iT[0][0]*cb*cl + iT[0][1]*cb*sl + iT[0][2]*sb;
  yeq = iT[1][0]*cb*cl + iT[1][1]*cb*sl + iT[1][2]*sb;
  sd = iT[2][0]*cb*cl + iT[2][1]*cb*sl + iT[2][2]*sb;
  cd = hypot(xeq, yeq);
  (*ad)[1] = atan2(sd,cd);
  //(*ad)[1] = asin(sd);
  //cd = cos((*ad)[1]);

  ca = xeq/cd;
  sa = yeq/cd;
  (*ad)[0] = atan2(sa,ca);
  if((*ad)[0]<0){(*ad)[0]+=twoPi;}

  (*ad)[2]=ca; (*ad)[3]=sa;   (*ad)[4]=cd; (*ad)[5]=sd;
}

void coords::ad2ecl(double a, double d, double* lambda, double* beta)
{
  vector<double> lambet;
  ad2ecl(a,d,&lambet);
  *lambda = lambet[0]; *beta = lambet[1];
}

void coords::ad2ecl(double a, double d, vector<double>* lambet)
{
  double cl, sl, cb, sb, cd, sd, ca, sa;
  lambet->resize(6);

  ca=cos(a); sa=sin(a);
  cd=cos(d); sd=sin(d);

  sb = T_ad_ecl[2][0]*cd*ca + T_ad_ecl[2][1]*cd*sa + T_ad_ecl[2][2]*sd;
  cb = sqrt(1.0 - sb*sb);
  (*lambet)[1] = atan2(sb,cb);
  //cb = cos((*lambet)[1]);

  cl = (T_ad_ecl[0][0]*cd*ca + T_ad_ecl[0][1]*cd*sa + T_ad_ecl[0][2]*sd)/cb;
  sl = (T_ad_ecl[1][0]*cd*ca + T_ad_ecl[1][1]*cd*sa + T_ad_ecl[1][2]*sd)/cb;
  (*lambet)[0] = atan2(sl,cl);
  if((*lambet)[0]<0){(*lambet)[0]+=twoPi;}
  if((*lambet)[0]>twoPi){(*lambet)[0]-=twoPi;}

  (*lambet)[2]=cl; (*lambet)[3]=sl;   (*lambet)[4]=cb; (*lambet)[5]=sb;
}

void coords::ecl2ad(double lambda, double beta, double* a, double* d)
{
  vector<double> ad;
  ecl2ad(lambda,beta,&ad);
  *a = ad[0]; *d = ad[1];
}

void coords::ecl2ad(double lambda, double beta, vector<double>* ad)
{
  double cl, sl, cb, sb, cd, sd, ca, sa;
  ad->resize(6);

  cl=cos(lambda); sl=sin(lambda);
  cb=cos(beta); sb=sin(beta);

  sd = iT_ad_ecl[2][0]*cb*cl + iT_ad_ecl[2][1]*cb*sl + iT_ad_ecl[2][2]*sb;
  cd = sqrt(1.0-sd*sd);
  (*ad)[1] = atan2(sd,cd);
  //(*ad)[1] = asin(sd);
  //cd = cos((*ad)[1]);

  ca = (iT_ad_ecl[0][0]*cb*cl + iT_ad_ecl[0][1]*cb*sl + iT_ad_ecl[0][2]*sb)/cd;
  sa = (iT_ad_ecl[1][0]*cb*cl + iT_ad_ecl[1][1]*cb*sl + iT_ad_ecl[1][2]*sb)/cd;
  (*ad)[0] = atan2(sa,ca);
  if((*ad)[0]<0){(*ad)[0]+=twoPi;}
  if((*ad)[0]>twoPi){(*ad)[0]-=twoPi;}

  (*ad)[2]=ca; (*ad)[3]=sa;   (*ad)[4]=cd; (*ad)[5]=sd;
}



void coords::XYZ2Rlb(vector<double> XYZ,vector<double>* Rlb)
{
  Rlb->resize(3);
  (*Rlb)[0] = sqrt(sqr(XYZ[0])+sqr(XYZ[1])+sqr(XYZ[2]));
  (*Rlb)[1] = atan2(XYZ[1],XYZ[0]);
  (*Rlb)[2] = asin(XYZ[2]/(*Rlb)[0]);
}

void coords::XYZ2Rlb(double X, double Y, double Z, double* R, double* l, double* b)
{
  vector<double> XYZ(3), Rlb(3);
  XYZ[0] = X;
  XYZ[1] = Y;
  XYZ[2] = Z;
  XYZ2Rlb(XYZ,&Rlb);
  *R = Rlb[0];
  *l = Rlb[1];
  *b = Rlb[2];
}

void coords::Rlb2XYZ(vector<double> Rlb, vector<double>* XYZ)
{
  XYZ->resize(3);
  (*XYZ)[2] = Rlb[0]*sin(Rlb[2]);
  double Dcb = Rlb[0]*cos(Rlb[2]);
  (*XYZ)[1] = Dcb*sin(Rlb[1]);
  (*XYZ)[0] = Dcb*cos(Rlb[1]);
}

void coords::Rlb2XYZ(double R, double l, double b, double* X, double* Y, double* Z)
{
  vector<double> Rlb(3), XYZ(3);
  Rlb[0] = R;
  Rlb[1] = l;
  Rlb[2] = b;
  Rlb2XYZ(Rlb,&XYZ);
  *X = XYZ[0];
  *Y = XYZ[1];
  *Z = XYZ[2];
}

void coords::uvw2rvmuadlb(double l, double b, double dist, double U, double V, double W, vector<double>* rvmuadlb)
{

  vector<double> ad;
  lb2ad(l,b,&ad);
  double a=ad[0];
  double d=ad[1];
  double ca=ad[2]; double sa=ad[3];
  double cd=ad[4]; double sd=ad[5];

  double UVW[3]; //U +ve towards gc for l,b=0,0. V +ve towards +ve Y. W +ve towards +ve Z
  UVW[0]=U; UVW[1]=V; UVW[2]=W;
  double pllx = 1.0/dist; //in whatever units are equivalent for d (kpc->mas)

  rvmuadlb->resize(5);

  double A[3][3], iB[3][3];
  A[0][0]=ca*cd; A[0][1]=-sa; A[0][2]=-ca*sd;
  A[1][0]=sa*cd; A[1][1]=ca;  A[1][2]=-sa*sd;
  A[2][0]=sd;    A[2][1]=0;   A[2][2]=cd;

  for(int i=0;i<3;i++)
    {
      for(int j=0;j<3;j++)
		{
		  iB[j][i]=0;
		  for(int k=0;k<3;k++)
			{
			  //we want the inverse = transpose/determinant; detB=1
			  iB[j][i] += T[i][k]*A[k][j];
			}
		}
    }

  for(int i=0;i<3;i++)
    {
      (*rvmuadlb)[i]=0;
      for(int j=0;j<3;j++)
		{
		  (*rvmuadlb)[i]+=iB[i][j]*UVW[j];
		}
    }

  for(int i=1;i<3;i++)
    {
      (*rvmuadlb)[i] *= pllx/4.74057;
    }


  // Convert the equatorial tangent velocity using the same Galactic rotation.
  muad2lb(a, d, (*rvmuadlb)[1], (*rvmuadlb)[2],
          &(*rvmuadlb)[3], &(*rvmuadlb)[4]);
  



  /*  double a2 = a + (*rvmuadlb)[1]*arcsec2rad/cd;
  double d2 = d + (*rvmuadlb)[2]*arcsec2rad;

  double l2, b2;
  ad2lb(a2,d2,&l2,&b2);
  
  (*rvmuadlb)[4] = (b2-b)/arcsec2rad;
  (*rvmuadlb)[3] = (l2-l)/arcsec2rad*cos(b);*/

  //result:
  //0 = rv; 1=mua; 2=mud; 3=mul; 4=mub;
  //proper motions are in mas yr-1 if dist in kpc
  //proper motions are in arcsec yr-1 if dist in pc

}

void coords::muad2lb(double a, double d, double mua, double mud, double* mul, double* mub)
{
  // Derive the Galactic-pole angles from T so this tangent transform follows
  // the same rotation as the position conversion.
  double angp = atan2(T[2][1], T[2][0]);
  if(angp < 0.0) angp += twoPi;
  const double cdngp = hypot(T[2][0], T[2][1]);
  const double sdngp = T[2][2];
  double cd = cos(d), sd = sin(d);
  double C1 = sdngp*cd - cdngp*sd*cos(a - angp);
  double C2 = cdngp*sin(a - angp);
  double cb = qAdd(C1,C2);

  *mul = (C1*mua + C2*mud)/cb;
  *mub = (-C2*mua + C1*mud)/cb;
}

void coords::mulb2ad(double l, double b, double mul, double mub, double* mua, double* mud)
{
  // Derive the node angle and pole latitude from T, matching muad2lb's frame.
  const double l0 = atan2(-T[0][2], T[1][2]);
  const double cdngp = hypot(T[2][0], T[2][1]);
  const double sdngp = T[2][2];
  double cl=cos(l-l0), sl=sin(l-l0);
  double cb=cos(b), sb=sin(b);
  double C1 = cb*sdngp - sb*cdngp*sl;
  double C2 = -cdngp*cl;
  double cd = qAdd(C1,C2);

  *mua = (C1*mul + C2*mub)/cd;
  *mud = (-C2*mul + C1*mub)/cd;
}

void coords::muad2ecl(double a, double d, double mua, double mud, double* mulambda, double* mubeta)
{
  //convert ad prop motion to lb prop motion (Based on refs within
  //Poleski (2013)

  //double ca=cos(a); double sa=sin(a);
  double cd=cos(d); double sd=sin(d);

  const double anep = 4.71238898038469; //270 deg
  const double cdnep = 0.39777696911260624;
  const double sdnep = 0.9174821430652418;
  
  double C1 = sdnep*cd - cdnep*sd*cos(a - anep);
  double C2 = cdnep*sin(a - anep);
  double cb = qAdd(C1,C2);

  *mulambda = (C1*mua + C2*mud)/cb;
  *mubeta = (-C2*mua + C1*mud)/cb;

}

void coords::muecl2ad(double a, double d, double mulambda, double mubeta, double* mua, double* mud)
{
  //convert ad prop motion to ad prop motion (Based on refs within
  //Poleski (2013)

  double lambda, beta;
  ad2ecl(a,d,&lambda,&beta);

  double lambda0 = 0.0;
  double cl=cos(lambda-lambda0); double sl=sin(lambda-lambda0);
  double cb=cos(beta); double sb=sin(beta);

  const double cdnep = 0.39777696911260624;
  const double sdnep = 0.9174821430652418;
  
  double C1 = cb*sdnep - sb*cdnep*sl;
  double C2 = -cdnep*cl;
  double cd = qAdd(C1,C2);

  *mua = (C1*mulambda + C2*mubeta)/cd;
  *mud = (-C2*mulambda + C1*mubeta)/cd;

}

double coords::fold(double x, double min, double max)
{
  double y;
  double tmp, diff, r;
  int N;

  if(max<min)
    {
      tmp=min;
      min=max;
      max=tmp;
    }

  if(min==max) return x; //wtf are you doing?

  if(x>=min && x<max) return x;

  r = max-min;

  if(x<min)
    {
      diff = x - min;
      N = floor(diff/r);
      y = x - N*r;
    }
  
  if(x>=max)
    {
      if(x==max) y=min;
      else
	{
	  diff = x - max;
	  N = ceil(diff/r);
	  y = x - N*r;
	}
    }

  return y;
}
