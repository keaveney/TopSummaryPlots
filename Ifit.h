#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include "TFile.h"
#include "TH2D.h"
#include "TVectorF.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TH2F.h"
#include "TMatrixD.h"

#include <iostream>
#include <fstream>

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>     /* atof */
#include "tdrstyle.C"
#include "CMS_lumi.C"


TGraphAsymmErrors * gr;
TGraphAsymmErrors * gr_cl;
TGraphAsymmErrors* g;

vector<string> yoda_files = {"TOP_16_011_DATA_NNLO.yoda","TOP_16_008_DATA_NNLO.yoda"};
vector<TGraphAsymmErrors*> graphs;


const int nbins = 15;


Double_t x[nbins],y[nbins],errory[nbins];
TVectorF vx(nbins), vy(nbins), vexl(nbins), vexh(nbins), veyl(nbins), veyh(nbins), clhi(nbins), cllo(nbins), line_y(nbins);

double outpar[2], err[2];



TFile * ff = new TFile("cov_file.root");
TH2F *cov = (TH2F*)ff->Get("inv_cov");


//options
bool custom_fit = true;
bool schmitt_fit = true;



Double_t func(float x,Double_t *par)
{
    // Double_t value=( (par[0]*par[0])/(x*x)-1)/ ( par[1]+par[2]*y-par[3]*y*y);
    
    Double_t value= par[0] + par[1]*x;
    return value;
}


Double_t dummy_func(float x, double par1, double par2 )
{
    // Double_t value=( (par[0]*par[0])/(x*x)-1)/ ( par[1]+par[2]*y-par[3]*y*y);
    
   // Double_t value= par[0] + par[1]*x;
    Double_t value= par1 + par2*x;
    return value;
}




TGraphAsymmErrors * read_YODA(std::string yoda_file){
    std::cout <<"in read_YODA   "<<   yoda_file <<std::endl;
    
    //vector<float> *vx, *vy, *vexl, *vexh, *veyl, *veyh;

    vector<string> tokens;
    string line;
    ifstream myfile (yoda_file);
    int row = 0;
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            tokens.clear();
            istringstream iss(line);
            copy(istream_iterator<string>(iss),
                 istream_iterator<string>(),
                 back_inserter(tokens));
            if (tokens.size() == 6){
                vx[row] = stof(tokens[0]);
            //    vexl[row] = stof(tokens[1]);
            //    vexh[row] = stof(tokens[2]);
                
                
                vexl[row] = 0.0;
                vexh[row] = 0.0;
                
                vy[row] = stof(tokens[3]);
                
                //veyl[row] = stof(tokens[3])*0.02;
                //veyh[row] = stof(tokens[3])*0.02;
                
               veyl[row] = stof(tokens[4]);
                veyh[row] = stof(tokens[5]);

                cout <<"VX = "<<  vx[row]  <<" VY "<<  vy[row] <<  " VEYL  " << veyl[row]  << " VEYH  " << veyh[row]<<endl;

                
                row++;
            }
        }
        myfile.close();

    }
    
    else cout << "Unable to open file";
    
    for (int i = 0; i < nbins; i++) cout <<" X val =  "<< vx[i];
    
    gr = new TGraphAsymmErrors (vx, vy, vexl, vexh, veyl, veyh);

    gr->Write("Ratio");
    return gr;
}



void create_cov_matrix(){


    TFile * cov_file =  new TFile("cov_file.root", "RECREATE");
    TH2F * cov = new TH2F("cov","cov", nbins,0.0,100,nbins,0,100);
    TH2F * inv_cov = new TH2F("inv_cov","inv_cov", nbins,0.0,100,nbins,0,100);
    TMatrixD mat =  TMatrixD(nbins,nbins);
    
    double covariance;
    double off_diag1;
    double off_diag2;
    double off_diag3;
    
    if(schmitt_fit == true){
         off_diag1 = 0.5;
         off_diag2 = 0.4;
         off_diag3 = 0.3;
    }else{
        off_diag1 = 0.0;
        off_diag2 = 0.0;
        off_diag3 = 0.0;
    }
    
    for (int i = 0; i < nbins; i++){
        for (int j = 0; j < nbins; j++){
            //covariance = (vy[i] - veyl[i]) * (vy[j] - veyl[j]);
            covariance = (veyl[i]) * (veyl[j]);
            
            cout <<" covanriance  = "<<veyl[i] <<  endl;
            
            
            if (i ==j){
               mat[i][j] = 1.0 * covariance;
               cov->SetBinContent(i+1,j+1,1.0*covariance);
            }
            else if (fabs(i - j) == 1.0) {
             cov->SetBinContent(i+1,j+1, off_diag1*covariance);
             mat[i][j] = off_diag1 * covariance;
            }
            else if (fabs(i - j) == 2.0) {
                cov->SetBinContent(i+1,j+1,off_diag2*covariance);
                mat[i][j] = off_diag2 * covariance;
            }
            else if (fabs(i - j) == 3.0){
                cov->SetBinContent(i+1,j+1,off_diag3*covariance);
                mat[i][j] = off_diag3 * covariance;
                     }
            else{
                cov->SetBinContent(i+1,j+1,0.0);
                mat[i][j] = 0.0;
                     }
        }
    }
    
    Double_t det2;
    mat.Invert(&det2);
   // TMatrixD inv_mat(kInvert,mat);
    
    
    for (int i = 0; i < nbins; i++){
        for ( int j = 0; j < nbins; j++){
                inv_cov->SetBinContent(i+1,j+1,mat[i][j]);
        }
    }


    TCanvas * c4 = new TCanvas();
    cov->Draw("COLZ");
    c4->SaveAs("covariance_matrix.png");
    

    inv_cov->Write();
    cov->Write();
    cov_file->Write();
    cov_file->Close();


}






double dummy_fcn(Double_t intercept, Double_t slope)
{
    
    Int_t i, j;
    
    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta[nbins];


    for (i=0;i<nbins; i++) {
        delta[i]  = (y[i]-dummy_func(x[i],intercept, slope));
    }
    

    
    for (i=0;i<nbins; i++) {
        for (j=0;j<nbins; j++) {
            double corr_coff = cov->GetBinContent(i+1,j+1);
            if (schmitt_fit) {
                chisq += delta[i]*delta[j]*corr_coff;
            }
            else{
                if(i==j) {
                    double res  = (y[i]-dummy_func(x[i],intercept, slope));
                    res = res*res;
                    res = res/errory[i];
                    chisq += res;
                }
            }
        }
    }
    
    return chisq;
}


void draw_contour(int nPoints){
    cout <<" drawing contour  "<< vx[0] <<endl;
    int scanp = nPoints*nPoints;
    
    TGraph2D *g2 = new TGraph2D(scanp);
    TGraph2D *g2_cont = new TGraph2D(scanp);
    
    vector<pair<float, float>> cl_band_vals;
    std::pair <float,float> vals;

    
    cout <<"out pars  HHH = = = "<<    outpar[0]    << "  "<< err[0] <<"  "<<  outpar[1]    << "  "<< err[1]  <<endl;
    double intercept = 0.0;
    double intercept_lo = (outpar[0] - (3*err[0]));
    double intercept_hi = (outpar[0] + (3*err[0]));
    
    double slope = 0.0;
    double slope_lo = (outpar[1] - (3*err[1]));
    double slope_hi = (outpar[1] + (3*err[1]));

    int i = 0;
    double min_chi2 = 1000000.00;
    double running_chi2 =0.0;
    
    for (intercept = intercept_lo; intercept < intercept_hi; intercept = intercept+(err[0]/100.0)){
        for (float slope = slope_lo; slope < slope_hi; slope = slope+(err[1]/100.0)){
            running_chi2  = dummy_fcn(intercept,slope);
            g2->SetPoint(i, intercept, slope, running_chi2);
            i++;
            if (running_chi2 < min_chi2 ) min_chi2 = running_chi2;
        }
    }

    cout <<" MIN CHI2  =   "<<   min_chi2<< endl;

    
    i = 0;
    for (intercept = intercept_lo; intercept < intercept_hi; intercept = intercept+(err[0]/100.0)){
        for (float slope = slope_lo; slope < slope_hi; slope = slope+(err[1]/100.0)){
            
            running_chi2  = dummy_fcn(intercept,slope);

            if ((running_chi2) < ((min_chi2) + 1.0 )){
                
              //  cout <<" IN BAND!   "<< endl;
              //  cout <<" MIN CHI2  =   "<<   min_chi2<< " running chi2  = = "<<  running_chi2  << " intercept =  "<< intercept<< " slope " <<slope<<endl;
                
                g2_cont->SetPoint(i, intercept, slope, 1.0);
                
                vals = std::make_pair(intercept,slope);
                cl_band_vals.push_back(vals);
            } else{
                g2_cont->SetPoint(i, intercept, slope, 0.0);
            }
            i++;
        }
    }
    
    
    cout <<"CL band vals size "<<  cl_band_vals.size() << endl;
    
    double running_yval = -999.0;
    double max_yval = -999.0;
    
    for (i = 0; i < nbins; i++){
        line_y[i] =  dummy_func(vx[i], outpar[0], outpar[1]);
    }
    
    //loop over graph points and estimate edges of CL band
    for (int i = 0; i < nbins; i++){
    
        max_yval = -999.0;
        cout <<" "<<endl;
        cout <<" Point  "<< i <<endl;
        for(int val = 0; val < cl_band_vals.size(); val++){
            

        
            running_yval = dummy_func(vx[i], cl_band_vals[val].first, cl_band_vals[val].second);
            

          //  cout <<"CL band vals "<<  cl_band_vals[val].first << "   "<<  cl_band_vals[val].second << endl;
          //  cout <<"Running yval = = = "<<  running_yval  << endl;
            
            if (running_yval > max_yval) max_yval = running_yval;
            
        }
        
        clhi[i] = fabs(max_yval - line_y[i]) ;
        cllo[i] = fabs(max_yval - line_y[i]);
    
    //    cout<<"   "<<endl;
       cout <<"Max yval = = = "<<  max_yval  << endl;
      cout<<"   "<<endl;
        
    }
    
    TCanvas * c4 = new TCanvas();

    gr_cl = new TGraphAsymmErrors (vx, line_y, vexl, vexh, cllo, clhi);

    
    gr_cl->Draw("ALE3");
    g->Draw("PSAME");
    gr_cl->SetFillColor(kRed);
    gr_cl->SetFillStyle(3002);
    gr_cl->GetXaxis()->SetRangeUser(0.0,800.0);
    gr_cl->GetYaxis()->SetRangeUser(0.6,1.4);
    gr_cl->GetHistogram()->GetXaxis()->SetTitle("Top p_{T}");
    gr_cl->GetHistogram()->GetYaxis()->SetTitle("DATA/NNLO");
    
    c4->SaveAs("Gr_CL.png");


    
    TCanvas * c3 = new TCanvas();
    g2->GetHistogram()->GetXaxis()->SetTitle("Intercept");
    g2->GetHistogram()->GetYaxis()->SetTitle("Slope");
    g2->Draw("COLZ");
    g2_cont->Draw("BOXSAME");
    g2_cont->SetMarkerStyle(2);
    c3->Write();
    c3->SaveAs("contour.png");

    TFile * fc = new TFile("chi2.root","RECREATE");
    g2->Write();
    fc->Write();
    
}



