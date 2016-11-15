//      f_nom->Draw("SAME");
//    f_nom->Write("fit_nocorrs");
//   g->Write("graph");
    


// g->SetMarkerStyle(5);
// g->SetMarkerSize(0.7);
// g->Draw("psame");
    
    
    
    
    
    
/*
  TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
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
  gMinuit->mnprin(3,amin);
    
//  gMinuit->mnexcm("MINOS",arglist,1,ierflg) ;
    
  double outpar[2], err[2], outpar_A[2], outpar_B[2];
    
  for (int i=0; i<2; i++){
       gMinuit->GetParameter(i,outpar[i],err[i]);
      
      if (i == 0){
          outpar_A[i] = outpar[i] + err[i];
          outpar_B[i] = outpar[i] - err[i];
      } else if(i == 1) {
          outpar_A[i] = outpar[i] - err[i];
          outpar_B[i] = outpar[i] + err[i];
      }
      }
    
    
  TGraph *gr0 = (TGraph *)gMinuit->Contour(80,0,1);
  TGraph *gr_fit = (TGraph*)gMinuit->GetPlot();

  TCanvas * c =  new TCanvas();

  TF1 *f_nom = new TF1("f_nom","[0] + [1]*x",0.0,800.0);
  f_nom->SetParameters(outpar);
    
  TF1 *f_up = new TF1("f_up","[0] + [1]*x",0.0,800.0);
  f_up->SetParameters(outpar_A);
    
  TF1 *f_down = new TF1("f_down","[0] + [1]*x",0.0,800.0);
  f_down->SetParameters(outpar_B);
    
    
    
    
    
    
  g->Draw("APE1");
  f_nom->Draw("SAME");
  f_up->Draw("SAME");
  f_down->Draw("SAME");
    
  gr_fit->Write();
  g->Write();
  gr0->Write();
  f->Write();
 
     
    
    TF1 *f_nom = new TF1("f_nom","[0] + [1]*x",0.0,800.0);
    f_nom->SetParameters(outpar);
*/
