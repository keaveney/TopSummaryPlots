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


#include "Ifit.h"

#include "TMinuit.h"
#include "TBackCompFitter.h"


using namespace std;


//______________________________________________________________________________

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Int_t i, j;
    
  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta[nbins];
    
    for (i=0;i<nbins; i++){
        delta[i]  = (y[i]-func(x[i],par));
    }
    
      for (i=0;i<nbins; i++) {
          for (j=0;j<nbins; j++) {
              double corr_coff = cov->GetBinContent(i+1,j+1);
              if (schmitt_fit) {
               chisq += delta[i]*delta[j]*corr_coff;
              }
              else {
                  if(i == j){
                      double res  = (y[i]-func(x[i],par));
                      res = res*res;
                      res = res/errory[i];
                      chisq += res;
                  }
                  }
     }
      }
    
  f = chisq;
}


//______________________________________________________________________________
int Ifit()
{
    setTDRStyle();

    
  TFile * f = new TFile("fit.root", "RECREATE");
    
    g = read_YODA("COMBINED_DATA_NNLO.yoda");
    
    
  // extract info from graph and fill arrays
  for (int i =0; i < nbins; i++){
    g->GetPoint(i,x[i],y[i]);
    errory[i] = g->GetErrorYhigh(i);
    }
    
  create_cov_matrix();
    
  TVirtualFitter * minuit = TVirtualFitter::Fitter(0,2);
  minuit->SetFCN(fcn);

  double arglist[10];
  arglist[0] = 1;
  double chi2, edm, errdef;
  int nvpar, nparx;
    
  // minimize
  arglist[0] = 100; // number of function calls
  arglist[1] = 1.0; // tolerance

    
  TF1 *f_nom = new TF1("f_nom","[0] + [1]*x",0.0,800.0);
    
  if(custom_fit){
      
       minuit->SetParameter(0, "Intercept",    0,     0.1, 0,0);
       minuit->SetParameter(1, "Slope",    0,     0.1, 0,0);
       minuit->ExecuteCommand("SET PRINT", arglist, 1);
       minuit->ExecuteCommand("MIGRAD", arglist, 0);
       minuit->ExecuteCommand("MINOS", arglist, 0);
       minuit->PrintResults(1,1.0);
       minuit->GetStats(chi2,edm,errdef,nvpar,nparx);

      
      for (int i=0; i<2; i++){
          gMinuit->GetParameter(i,outpar[i],err[i]);
      }
      
      cout <<"out pars  CCC = = = "<<    outpar[0]    << "  "<< err[0] <<"  "<<  outpar[1]    << "  "<< err[1]  <<endl;

   f_nom->SetParameters(outpar);
   g->SetMarkerStyle(23);
   g->Draw("APE1");
   f_nom->Draw("same");
      
   //Create a TGraphErrors to hold the confidence intervals
  TGraphErrors *grint = new TGraphErrors(9);
  grint->SetTitle("Fitted line with .95 conf. band");
      
  for (int i=0; i<9; i++) {
        grint->SetPoint(i, g->GetX()[i], 0);
      }
      

  grint->SetFillStyle(3001);
  grint->SetFillColor(kRed);
  grint->SetLineColor(kRed);
  grint->Draw("E3same");
    
  }
    else{

        g->Fit(f_nom);
        g->Draw("APE1");

        
  //Create a TGraphErrors to hold the confidence intervals
  TGraphErrors *grint = new TGraphErrors(9);
  grint->SetTitle("Fitted line with .95 conf. band");
    
  for (int i=0; i<9; i++) {
      grint->SetPoint(i, g->GetX()[i], 0);
    }
    
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint);
   grint->SetFillStyle(3001);
   grint->SetFillColor(kRed);
   grint->SetLineColor(kRed);

   grint->Draw("E3same");
    }
    
    
    draw_contour(100);

    
   cout <<"CHI2 = "<<  chi2   << endl;

  return 1;

}





