#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH1.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TFile.h"


void ConfidenceIntervals()
{
  //Illustrates TVirtualFitter::GetConfidenceIntervals
  //This method computes confidence intervals for the fitted function
  //Author: Anna Kreshuk

  TCanvas *myc = new TCanvas("myc",
			     "Confidence intervals on the fitted function",1200, 500);
 // myc->Divide(3,1);

  /////1. A graph
  //Create and fill a graph
  Int_t ngr = 100;
    
TFile * ff = new TFile("fit.root");

  TGraphErrors *gr = (TGraphErrors*)ff->Get("Ratio");
    
  //  for (int i=0;i<9; i++) {
  //      gr->SetPointError(i,0.00001,  gr->GetErrorY(i));
  //  }
    
    
//TH2F *cov = (TH2F*)ff->Get("inv_cov");

//  gr->SetName("GraphNoError");
  Double_t x, y;
  Int_t i;

    
//  for (i=0; i<ngr; i++){
//    x = gRandom->Uniform(-1, 1);
//    y = -1 + 2*x + gRandom->Gaus(0, 1);
//    gr->SetPoint(i, x, y);
//  }
    
    
    
  //Create the fitting function
  TF1 *fpol = new TF1("fpol", "pol1", -1, 1);
  fpol->SetLineWidth(2);
  gr->Fit(fpol, "Q");

  //Create a TGraphErrors to hold the confidence intervals
  TGraphErrors *grint = new TGraphErrors(ngr);
  grint->SetTitle("Fitted line with .95 conf. band");
  for (i=0; i<ngr; i++)
    grint->SetPoint(i, gr->GetX()[i], 0);
  //Compute the confidence intervals at the x points of the created graph
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint);
  //Now the "grint" graph contains function values as its y-coordinates
  //and confidence intervals as the errors on these coordinates
  //Draw the graph, the function and the confidence intervals
  grint->SetFillColor(kRed);
  grint->Draw("ae3");
  gr->SetMarkerStyle(5);
  gr->SetMarkerSize(0.7);
  gr->Draw("psame");

  
  }
