//
//   Example of a program to fit non-equidistant data points
//   =======================================================
//
//   The fitting function fcn is a simple chisquare function
//   The data consists of 5 data points (arrays x,y,z) + the errors in errorsz
//   More details on the various functions or parameters for these functions
//   can be obtained in an interactive ROOT session with:
//    Root > TMinuit *minuit = new TMinuit(10);
//    Root > minuit->mnhelp("*")  to see the list of possible keywords
//    Root > minuit->mnhelp("SET") explains most parameters
//Author: Rene Brun

#include "TMinuit.h"

Float_t x[9],y[9],errory[9];

//______________________________________________________________________________
Double_t func(float x,Double_t *par)
{
 // Double_t value=( (par[0]*par[0])/(x*x)-1)/ ( par[1]+par[2]*y-par[3]*y*y);
    
  Double_t value= par[0] + par[1]*x;
  return value;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  const Int_t nbins = 9;
  Int_t i;

  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta;
  for (i=0;i<nbins; i++) {
 //   delta  = (z[i]-func(x[i],y[i],par))/errorz[i];
      delta  = (y[i]-func(x[i],par))/errory[i];
      chisq += delta*delta;
  }
  f = chisq;
}

//______________________________________________________________________________
void Ifit()
{

  // the x values
  x[0]=22.5;
  x[1]=67.5;
  x[2]=112.5;
  x[3]=157.5;
  x[4]=202.5;
  x[5]=247.5;
  x[6]=292.5;
  x[7]=357.5;
  x[8]=600.5;
    
  // the y values
  y[0]=1.0/0.986;
  y[1]=1.0/0.974;
  y[2]=1.0/0.0992;
  y[3]=1.0/1.018;
  y[4]=1.0/1.042;
  y[5]=1.0/1.062;
  y[6]=1.0/1.023;
  y[7]=1.0/1.036;
  y[8]=1.0/1.083;

    
// The errors on y values
  Float_t error = 0.01;
  errory[0]=0.0842;
  errory[1]=0.0497;
  errory[2]=0.0228;
  errory[3]=0.0584;
  errory[4]=0.0638;
  errory[5]=0.0651;
  errory[6]=0.114;
  errory[7]=0.162;
  errory[8]=0.146;
    

  TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 5 params
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  // Set starting values and step sizes for parameters
  static Double_t vstart[2] = {0.0,0.0};
  static Double_t step[2] = {0.001, 0.001};
  gMinuit->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
//  gMinuit->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);
//  gMinuit->mnparm(3, "a4", vstart[3], step[3], 0,0,ierflg);

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);

}
