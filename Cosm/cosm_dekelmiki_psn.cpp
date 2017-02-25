/*****************
; This program is a cpp translation of the code written by Dekel .
; Function that returns the cumulative number density, the characteristic mass
; and sigma(M,z) using the ST approximation. It follows the 
; the formalism summarized in Dekel&Birnboim 2006 appendix.
;
;******************************/

//include standard libraries
#include <sstream>
#include <iostream> 
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string.h> 
#include <fstream>
#include <vector>
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <algorithm>

using namespace std;
//include header
double ez2();
void signor();
double sig(double rs);
double su(double rk, double rs);
double pk(double rk);
double w(double x);
double ommz();
double omlz();
double dz();
double delta();
double sigm(double em);
double emstar();

//Global variables
double omm,oml,h,sigma8,signorm,rho0,z,pi;


/*******************************
Start with a few subroutine
********************************/

//ez2: checked
double ez2(){
  return oml+(1.-oml-omm)*(1.+z)*(1.+z)+omm*(1.+z)*(1.+z)*(1.+z);
  
}


//Compute the approximated power spectrum: checked
double pk(double rk){
  double q=rk/(omm*h*h); 
  double sog=1.+3.89*q+(16.1*q)*(16.1*q)+pow((5.46*q),3.)+pow((6.71*q),4.);
  double t=log(1.+2.34*q)/(2.34*q)/pow(sog,0.25);
  return signorm*signorm*rk*t*t;
}


//Compute the filter function: checked
double w(double x){
  return 3.*(sin(x)-x*cos(x))/pow(x,3.);
}

//Define the integrand of sigma(R): checked
double su(double rk, double rs){
  return rk*rk*pk(rk)*w(rk*rs)*w(rk*rs);
}


//compute sigma(R) integrating over power spectrum: checked
double sig(double rs){
  //set the boundary and precision
  double rklmin=-3.; 
  double rklmax=2.-log10(rs);
  double dkl=0.001;
  double rkmin=pow(10.,rklmin);
  double rkmax=pow(10.,rklmax);

  //make the integral with trapezoid rule
  double sum=0.5*dkl*(rklmin*su(rkmin,rs)+rklmax*su(rkmax,rs));
  double rk;
  for(double rkl=rklmin+dkl; rkl<rklmax; rkl=rkl+dkl){
    rk=pow(10.,rkl);
    sum=sum+dkl*rk*su(rk,rs);
  }

  //set final values
  double sig2=sum/2./(pi*pi)*log(10.);
  return sqrt(sig2);
}

//normalize to sigma8: checked
void  signor(){
  //call once before emstar
  double rs=8./h;
  signorm=1.;
  double sig8=sig(rs);
  signorm=sigma8/sig8;
  return;
}


//Cosmological constants with redhsift: checked
double ommz(){
  return omm*pow((1.+z),3.)/ez2();
}
double omlz(){
  return oml/ez2();
}


//Compute D(a): checked
double dz(){
  double omz=ommz();
  double olz=omlz();
  double gz=2.5*omz/(pow(omz,(4./7.))-olz+(1.+omz/2.)*(1.+olz/70.));
  double g0=2.5*omm/(pow(omm,(4./7.))-oml+(1.+omm/2.)*(1.+oml/70.));
  return gz/(g0*(1.+z));
}


//Delta z: checked
double delta(){
  double a=1./(1.+z);
  double ommz=omm/(a*a*a)/ez2();
  double x=ommz-1.;
  return (18.*pi*pi+82.*x-39.*x*x)/(1.+x);
}

 
//Get sigma m: checked
double sigm(double em){
  double rs=pow((em/(4.*pi/3.)/rho0),(1./3.));
  return sig(rs);
}

//mstar: checked
double emstar(){    
  double un=1.;  // nu
  //double aa=1./(1.+z)/((delta()/200.*omm/0.3*(h/0.7)*(h/0.7)),(1./3.));  
  double sigc=1.686/dz()/un;
  double em0=130.*pow((un*dz()),8.)/100.; //low bound good z=0-6. For z=10 do **10(?) 
  double em=em0;
  double dm=em0;
  double rs,sigma;
  do{
    em=em+dm;
    rs=0.143*pow((200.*em/(omm/0.3*(h/0.7)*(h/0.7))),(1./3.));
    sigma=sig(rs);
  }while(sigma > sigc);
  if(em == em0+dm){
    cout<<"em0 too big!"<<endl;
    exit(-1);
  }
  return em-0.5*dm;
}

int main(int argc, char* argv[]){

  //parse the input
  if(argc < 8){
    cout<<"Missing parameters!!!";
    exit(-1);
  }

  //set common variable
  omm=(double)atof(argv[4]);
  z=(double)atof(argv[3]);
  h=(double)atof(argv[5]);
  sigma8=(double)atof(argv[6]);
  pi=4*atan(1.);
  
  //set mass boundary and precision in integral
  double eminl=(double)(atof(argv[1])-11);
  double emaxl=(double)(atof(argv[2])-11);
  double deml=(double)atof(argv[7]);
    
  //set other cosmology
  //double fb=0.1;
  double delc=1.686;
    
  //set other global
  oml=1.0-omm;
  rho0=2.76*omm*h*h;          //units of 10^11 M_solar / (Mpc)^3  comoving
  //double h0=0.0716*(h/0.7);   // 1/Gyr
  

  //normalize by sigma8 and update global
  signor(); //checked. It works


  //integrated to get the total number density
  double enn=0.;
  
  
  double em, un, fst, eml1, eml2, em1, em2, sigp, dundm, en, eml;    
  for(eml=eminl+deml;eml<=emaxl;eml=eml+deml){
    em=pow(10.,eml);
    un=delc/dz()/sigm(em);
    un=0.841*un;
    fst=0.322*(1.+1./pow(un,0.6))*0.841;
    eml1=eml-0.5*deml;
    eml2=eml+0.5*deml;
    em1=pow(10.,eml1);
    em2=pow(10.,eml2);
    sigp=(sigm(em2)-sigm(em1))/deml; // note d/dlogM 
    dundm=-delc/dz()*sigp/(sigm(em)*sigm(em));
    en=fst*sqrt(2./pi)*rho0/em*dundm*exp(-0.5*un*un); 
    enn=enn+en*deml;
  }
  
  
  //apply normalsation
  eml=eminl;
  em=pow(10.,eml);
  un=delc/dz()/sigm(em);
  un=0.841*un;
  fst=0.322*(1.+1./pow(un,0.6))*0.841;
  eml1=eml-0.5*deml;
  eml2=eml+0.5*deml;
  em1=pow(10.,eml1);
  em2=pow(10.,eml2);
  sigp=(sigm(em2)-sigm(em1))/deml;
  dundm=-delc/dz()*sigp/(sigm(em)*sigm(em));
  en=fst*sqrt(2./pi)*rho0/em*dundm*exp(-0.5*un*un);
  enn=enn+en*deml*0.5;
  
  //set the return so that spawn can parse it
  //number density, mstar
  cout<<enn<<"W"<<emstar()<<"Z"<<endl;
  exit(0);
}

