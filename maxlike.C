//  example of a simple max likelihood estimation
//  inspired in the examples from Cowan, chap.6 and in root
//  tutorial Ifit.C
//   - generate fake data according to a given function
//   - get the parameters with the max likelihood method
//   - plot histograms

//                  jtmn setp 2003
//
//  comments from Ifit.C:
//   More details on the various functions or parameters for these functions
//   can be obtained in an interactive ROOT session with:
//    Root > TMinuit *minuit = new TMinuit(10);
//    Root > minuit->mnhelp("*")  to see the list of possible keywords
//    Root > minuit->mnhelp("SET") explains most parameters
//

#include "TMinuit.h"
// global variables
const Int_t Npoints = 2000;
Double_t xp[Npoints];
Double_t xmin = -0.95;
Double_t xmax = +0.95;
Double_t alpha = 0.5; 
Double_t beta  = 0.5;
//______________________________________________________________________________
Double_t fgen(Double_t x)
{
  // function normalized in the interval -1,1

  Double_t value = (1.0 + alpha*x + beta*x*x)/(2.0 + 2.0*beta/3.0);
  return value;
}
//______________________________________________________________________________
Double_t func(Double_t x,Double_t *par)
{
  Double_t xmin2 = TMath::Power(xmin,2);
  Double_t xmin3 = TMath::Power(xmin,3);
  Double_t xmax2 = TMath::Power(xmax,2);
  Double_t xmax3 = TMath::Power(xmax,3);
  Double_t value = (1.0 + par[0]*x + par[1]*x*x)/ 
    ((xmax-xmin) + par[0]*(xmax2-xmin2)/2.0 + par[1]*(xmax3-xmin3)/3.0);
  return value;
}
//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // f is the function that MINUIT really minimizes (so define it with - in front
  // to maximize
  Int_t i;
  Double_t logl = 0;
  Double_t temp;
  for (i=0; i<Npoints;i++) {
    temp = TMath::Log(func(xp[i],par));
    logl += temp;
  }
  f = -logl;
}

//______________________________________________________________________________
void maxlike()
{
  TRandom *r3 = new TRandom3();
  r3->SetSeed(0);
  // generate the fake data using the Von Neumann's acceptance-rejection method
  Int_t naccept=0;
  while  (naccept < Npoints)
    { 
      // x goes from -1 to +1
      Double_t x = 2.0*(r3->Rndm())-1.0;
      // y goes from 0 to fgen(1) (for alpha,beta >0)
      Double_t y = fgen(1)*(r3->Rndm());
      if (y<=fgen(x)) {
	// accept the point if between xmax and xmin
        if ((x<xmax) && (x>xmin)){
	  xp[naccept]=x;
	  naccept++;
	}
      }
    }

  // show histogram
  TCanvas *fake = new TCanvas("fake"," Hist  ", 600,600);
  TH1F *h  = new TH1F("h"," ", 20,-1.0,1.0);
  for (Int_t i =0; i<Npoints; i++){
    h->Fill(xp[i]);
  }
  h->DrawNormalized();

  // prepare maximization (in fact minimization)
  TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;

  // for max likelihood = 0.5, for chisq = 1.0
  arglist[0] = 0.5;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  // Set starting values and step sizes for parameters
  static Double_t vstart[2] = {1.,1.,};
  static Double_t step[2] = {0.001 , 0.001};
  gMinuit->mnparm(0, "alpha", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "beta",  vstart[1], step[1], 0,0,ierflg);

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);
   
  //now make the contours (taken from fitcont.C)

  TCanvas *c2 = new TCanvas("c2","contours",10,10,600,600);
  //Get contour for parameter 0 versus parameter 1  for ERRDEF=2 
  gMinuit->SetErrorDef(2);
  TGraph *gr2 = (TGraph*)gMinuit->Contour(80,0,1);
  gr2->SetFillColor(42);
  gr2->SetTitle("Exemplo do Minuit -- max. verossimilhanca");
  gr2->Draw("alf");
  //Get contour for parameter 0 versus parameter 1 for ERRDEF=1  
  gMinuit->SetErrorDef(1);
  TGraph *gr1 = (TGraph*)gMinuit->Contour(80,0,1);
  gr1->SetFillColor(38);
  gr1->Draw("lf");
  // put the star as the TRUE value
  trueP = new TMarker(alpha,beta,8);
  trueP->SetMarkerColor(kRed);
  trueP->SetMarkerSize(2);
  trueP->Draw();
  // get the estimated values
  Double_t alpha_est, alpha_est_err;
  Double_t beta_est, beta_est_err;
  gMinuit->GetParameter(0,alpha_est,alpha_est_err);
  gMinuit->GetParameter(1,beta_est,beta_est_err);
  // put the square as the estimated value
  estP = new TMarker(alpha_est,beta_est,8);
  estP->SetMarkerColor(kBlack);
  estP->SetMarkerSize(2);
  estP->Draw();


}
