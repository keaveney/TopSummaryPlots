#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMathText.h"
#include "TText.h"
#include "TLine.h"
#include "TLatex.h"


#include "TFile.h"
#include "TH2D.h"
#include "TVectorF.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMatrixD.h"

#include <iostream>
#include <fstream>
#include <math.h>       /* exp */


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

TGraphAsymmErrors* g_temp;

TF1 *f_nom1;

Double_t res;
Double_t f;

vector<string> yoda_files = {"TOP_16_011_DATA_NNLO.yoda","TOP_16_008_DATA_NNLO.yoda"};
//vector<string> yoda_files = { "TOP_16_011_DATA_PWHGP8.yoda", "TOP_16_008_DATA_PWHGP8.yoda"};


vector<TGraphAsymmErrors> graphs;


//style options for TDRstyle
int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
int iPos = 11;


string fit = "exp";

const int nbins = 15;

//global array needed for objective fitting function and
//writing of covariance matrix
Double_t x[nbins],y[nbins],errory[nbins];

double outpar[2], err[2];


TFile * ff = new TFile("cov_file.root");
TH2F *cov = (TH2F*)ff->Get("inv_cov");


//options
bool custom_fit = true;
bool schmitt_fit = true;
bool BCC = true;


Double_t func(float x,Double_t *par)
{
    // Double_t value=( (par[0]*par[0])/(x*x)-1)/ ( par[1]+par[2]*y-par[3]*y*y);
    
    //Double_t value= par[0] + par[1]*x;
    Double_t value;
    
    if (fit == "linear"){
        
        value= par[0] + par[1]*x;
        
    }else if(fit == "exp"){
        
        value = exp(par[0] + par[1]*x);
        
    }
    return value;
}


Double_t dummy_func(float x, double par1, double par2 )
{
    // Double_t value=( (par[0]*par[0])/(x*x)-1)/ ( par[1]+par[2]*y-par[3]*y*y);
    
   // Double_t value= par[0] + par[1]*x;
//Double_t value= par1 + par2*x;
    
    Double_t value;

    if (fit == "linear"){
        
        value= par1 + par2*x;
        
    }else if(fit == "exp"){
        
        value = exp(par1 + par2*x);
        
    }
    
    

    return value;
}



TGraphAsymmErrors * read_YODA(std::string yoda_file){
    std::cout <<"in read_YODA   "<<   yoda_file <<std::endl;
    
    int nbins = 0;

    
    if ((yoda_file == "TOP_16_011_DATA_NNLO.yoda")   || (yoda_file == "TOP_16_011_DATA_PWHGP8.yoda")   ){
    nbins = 6;
        
    } else if((yoda_file == "TOP_16_008_DATA_NNLO.yoda")    || (yoda_file == "TOP_16_008_DATA_PWHGP8.yoda")   ){
        nbins = 9;
    }
    
    else if(yoda_file == "COMBINED_DATA_NNLO.yoda"){
        
        nbins = 15;
        
    }
    
    
    Double_t xsecs_ll[6] = {4.05e-03, 6.13e-03, 3.30e-03, 1.03e-03, 2.12e-04, 3.91e-05};
    Double_t xsecs_ljets[9] = {0.00296538, 0.00632151, 0.00558930, 0.00348763, 0.00185964, 0.00093842 , 0.00049746 , 0.000201128, 0.0000205860 };
    


    
    TVectorF vx(nbins), vy(nbins), vexl(nbins), vexh(nbins), veyl(nbins), veyh(nbins), clhi(nbins), cllo(nbins), line_y(nbins);
    Double_t bin_edges[nbins+1];
    
    
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
                vexl[row] = stof(tokens[1]);
                vexh[row] = stof(tokens[2]);
                bin_edges[row] = (stof(tokens[0]) - stof(tokens[1]));
                
             //   cout << "Bin " << row <<" , bin centre  "<< stof(tokens[0])  <<"bin error "<<  stof(tokens[1]) <<" bin edge " << bin_edges[row] << endl;
                
              //  vexl[row] = 0.0;
              //  vexh[row] = 0.0;
                
                vy[row] = stof(tokens[3]);
                
                //veyl[row] = stof(tokens[3])*0.02;
                //veyh[row] = stof(tokens[3])*0.02;
                
                veyl[row] = stof(tokens[4]);
                veyh[row] = stof(tokens[5]);

               // cout <<"VX = "<<  vx[row]  <<" VY "<<  vy[row] <<  " VEYL  " << veyl[row]  << " VEYH  " << veyh[row]<<endl;
                
                row++;
            }
        }
        myfile.close();
        bin_edges[nbins] = (vx[nbins-1] + vexh[nbins-1]);

        cout <<"closed file"<<endl;

    }
    
    else cout << "Unable to open file";
    
    //Calculate Bin-Centre Corrections (BCC)
    if (BCC){
        
        double raw_xsec;
        cout << "Calculating bin centre corrections...";
        
        TFile * f_fine = new TFile("emu_ttbarsignalplustau.root");
        TH1F  * h_fine = (TH1F*)f_fine->Get("VisGenToppT");

        double tt_xs = 1.0;
        double norm_factor = tt_xs/h_fine->Integral();

        h_fine->Scale(norm_factor);
        TCanvas* c_binning = new TCanvas();
        h_fine->Draw();
        

        for (int yval = 0; yval < nbins; yval++){
        
            double yval_res = 999.00, min_yval_res = 999.00;
            double corrected_xval = vx[yval];
            cout <<"  "<<endl;
            cout <<"BROAD BIN NUMBER   "<< yval <<endl;
            
            if ((yoda_file == "TOP_16_011_DATA_NNLO.yoda")    ||   (yoda_file == "TOP_16_011_DATA_PWHGP8.yoda")){
                raw_xsec = xsecs_ll[yval];
            } else if((yoda_file == "TOP_16_008_DATA_NNLO.yoda")  ||   (yoda_file == "TOP_16_008_DATA_PWHGP8.yoda")   ){
                raw_xsec = xsecs_ljets[yval];
            }
            
            TLine * l_bin = new TLine(0.0,raw_xsec,800.0,raw_xsec);
            l_bin->SetLineStyle(2);
            l_bin->SetLineColor(2);
            l_bin->Draw();
            
            double original_xval =vx[yval] ;
            
            for (int xbin_fine = 1; xbin_fine < h_fine->GetNbinsX(); xbin_fine++){
                if ((h_fine->GetBinCenter(xbin_fine)>(vx[yval]-vexl[yval]))&&(h_fine->GetBinCenter(xbin_fine)<(vx[yval]+vexh[yval]))){
                    
                yval_res = fabs(h_fine->GetBinContent(xbin_fine)  - raw_xsec);
                  // cout <<" looping:   FINE BIN NUMBER =  "<< xbin_fine << " , fine bin centre =   "<<  h_fine->GetBinCenter(xbin_fine) << " resdual  = "<<  yval_res << endl;

                if (yval_res < min_yval_res){
                    min_yval_res = yval_res;
                    corrected_xval = h_fine->GetBinCenter(xbin_fine);
                    vx[yval] = corrected_xval;
                    x[yval] = corrected_xval;
                }
                }
                
                
            }
            
            cout <<" min res = "<< min_yval_res  << " original X val  = "<<   original_xval  <<", corrected X val = "<<  corrected_xval <<endl;

        }
        }

    //Set X bin errors to 0 for convenience
    for (int xbin = 0; xbin < nbins; xbin++){
        vexl[xbin] = 0.0;
        vexh[xbin] = 0.0;
    }
    
    
    
    gr = new TGraphAsymmErrors(vx, vy, vexl, vexh, veyl, veyh);
    gr->SetName(yoda_file.c_str());
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
         off_diag1 = 0.3;
         off_diag2 = 0.2;
         off_diag3 = 0.0;
    }else{
        off_diag1 = 0.0;
        off_diag2 = 0.0;
        off_diag3 = 0.0;
    }
    
    for (int i = 0; i < nbins; i++){
        for (int j = 0; j < nbins; j++){
            //covariance = (vy[i] - veyl[i]) * (vy[j] - veyl[j]);
            covariance = (errory[i]) * (errory[j]);
            
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






double dummy_fcn(Double_t  intercept, Double_t slope)
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
    
    TVectorF vx(nbins), vy(nbins), vexl(nbins), vexh(nbins), veyl(nbins), veyh(nbins), clhi(nbins), cllo(nbins), line_y(nbins);

    
    int scanp = nPoints*nPoints;
    
    TGraph2D *g2 = new TGraph2D(scanp);
    TGraph2D *g2_cont = new TGraph2D(scanp);
    
    vector<pair<float, float>> cl_band_vals;
    std::pair <float,float> vals;
    
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
                
              // cout <<" IN BAND!   "<< endl;
               // cout <<" MIN CHI2  =   "<<   min_chi2<< " running chi2  = = "<<  running_chi2  << " intercept =  "<< intercept<< " slope " <<slope<<endl;
                
                g2_cont->SetPoint(i, intercept, slope, 1.0);
                
                vals = std::make_pair(intercept,slope);
                cl_band_vals.push_back(vals);
            } else{
                g2_cont->SetPoint(i, intercept, slope, 0.0);
            }
            i++;
        }
    }
    
    
    double running_yval = -999.0;
    double max_yval = -999.0;
    double even_x = 0.0;
    
    for (i = 0; i < nbins; i++){
        even_x = even_x + (i*(800.0/nbins));
        vx[i] = even_x;
        line_y[i] =  dummy_func(even_x, outpar[0], outpar[1]);

    }
    
    //loop over graph points and estimate edges of CL band
    for (int i = 0; i < nbins; i++){
    
        max_yval = -999.0;
        for(int val = 0; val < cl_band_vals.size(); val++){
            
            running_yval = dummy_func(vx[i], cl_band_vals[val].first, cl_band_vals[val].second);
            
            if (running_yval > max_yval) max_yval = running_yval;
            
        }
        
        clhi[i] = fabs(max_yval - line_y[i]) ;
        cllo[i] = fabs(max_yval - line_y[i]);

    }
    
    for(int i = 0 ; i < nbins; i++){
       // vx[i] = x[i];
        vy[i] = y[i];
        vexl[i] = 0.0;
        vexh[i] = 0.0;
        veyl[i] = errory[i];
        veyh[i] = errory[i];
            }

    
    gr_cl = new TGraphAsymmErrors (vx, line_y, vexl, vexh, cllo, clhi);

    
    TCanvas * c4 = new TCanvas();

    gr_cl->Draw("AE3");
    
    vector<int> colors = {2,4};
    vector<int> styles = {22,23};
    vector<double> sizes = {1.3,1.3};
    
    gr_cl->SetFillColor(kBlack);
    gr_cl->SetFillStyle(3002);
    gr_cl->GetXaxis()->SetRangeUser(0.0,690.0);
    gr_cl->GetYaxis()->SetRangeUser(0.5,1.5);
    gr_cl->GetHistogram()->GetXaxis()->SetTitle("Top p_{T} [GeV]");
    
    if (yoda_files[0].find("PWHGP8") != std::string::npos) {
        gr_cl->GetHistogram()->GetYaxis()->SetTitle("Data/Powheg + Pythia8");
    }
    else if(yoda_files[0].find("NNLO") != std::string::npos){
        gr_cl->GetHistogram()->GetYaxis()->SetTitle("Data/NNLO");
    }
    gr_cl->SetName("fit");

    if(fit== "exp"){
        f_nom1 = new TF1("f_nom1","exp([0] + (x*[1]))",0.0,800.0);
    
    } else if (fit == "linear"){
        f_nom1 = new TF1("f_nom1","[0] + (x*[1])",0.0,800.0);
    }
    
    f_nom1->SetParameters(outpar);
    f_nom1->SetLineColor(1);
    f_nom1->Draw("same");


    TLegend *leg = new TLegend(0.42,0.68,0.7,0.82);
    leg->SetTextFont(65);
    TLegendEntry* l1 = leg->AddEntry("TOP_16_011_NNLO_DATA.yoda" ,"dilepton (CMS-PAS-TOP-16-011)","E1p");
    l1->SetMarkerColor(2);
    l1->SetMarkerStyle(22);
    l1->SetLineColor(2);
    TLegendEntry* l2 = leg->AddEntry("TOP_16_008_NNLO_DATA.yoda","lepton + jets (arXiv:1610.04191, sub. to PRD)","E1p");
    l2->SetMarkerColor(4);
    l2->SetMarkerStyle(23);
    l2->SetLineColor(4);
    TLegendEntry* l3 = leg->AddEntry("fit","fit","lf");
    leg->SetBorderSize(0.0);
    leg->Draw();
   

    
    for (int i = 0 ; i < graphs.size(); i++){
        graphs[i].SetMarkerColor(colors[i]);
        graphs[i].SetMarkerStyle(styles[i]);
        graphs[i].SetMarkerSize(sizes[i]);
        graphs[i].SetLineColor(colors[i]);
        graphs[i].Draw("PSAME");
    }
    
    
   /*
    TLegend *leg = new TLegend(0.42,0.68,0.7,0.82);
    leg->SetTextFont(65);
    leg->AddEntry("TOP_16_011_NNLO_DATA.yoda" ,"dilepton (CMS-PAS-TOP-16-011)","E1p");
    leg->AddEntry("TOP_16_008_NNLO_DATA.yoda","lepton + jets (arXiv:1610.04191, sub. to PRD)","E1p");
    leg->AddEntry("fit","fit","lf");
    leg->SetBorderSize(0.0);
    leg->Draw();
    */
    
    std::string s = "ratio = e^{";
    
    std::ostringstream ss;
    ss <<  std::setprecision(4) << std::fixed <<outpar[0];
    std::string outpar1_str(ss.str());
    
    ss.str(std::string());
    double outpar2_s = outpar[1]*(-1.0);
    
    ss << std::setprecision(4) << std::fixed << outpar2_s;
    std::string outpar2_str(ss.str());

    std::string result_string = s + outpar1_str +" - " +outpar2_str + "  #times  p_{T}}";
    
    
    TText * leg_text = new TText();
    leg_text->SetTextAlign(23);
    leg_text->SetTextSize(0.045);
    leg_text->DrawText(400.0, 1.4, "Parton level");
    
    TLatex * l = new TLatex();
    l->SetTextAlign(23);
    l->SetTextSize(0.045);
    l->DrawLatex(250.0, 0.6, result_string.c_str());

    CMS_lumi( c4, iPeriod, iPos );

    c4->SaveAs("Gr_CL.png");

    TCanvas * c3 = new TCanvas();
    g2->GetHistogram()->GetXaxis()->SetTitle("P[0]");
    g2->GetHistogram()->GetYaxis()->SetTitle("P[1]");
    g2->Draw("COLZ");
    g2_cont->Draw("BOXSAME");
    g2_cont->SetMarkerStyle(2);
    
    double x1_intercept_lo,  x2_intercept_lo, y1_intercept_lo,  y2_intercept_lo;
    double x1_intercept_hi,  x2_intercept_hi, y1_intercept_hi,  y2_intercept_hi;
    double x1_slope_lo,  x2_slope_lo, y1_slope_lo,  y2_slope_lo;
    double x1_slope_hi,  x2_slope_hi, y1_slope_hi,  y2_slope_hi;

    
    x1_intercept_lo = outpar[0] - (err[0]);
    x2_intercept_lo = outpar[0] - (err[0]);
    y1_intercept_lo = outpar[1] - (err[1]);
    y2_intercept_lo = outpar[1] + (err[1]);
    
    x1_intercept_hi = outpar[0] + (err[0]);
    x2_intercept_hi = outpar[0] + (err[0]);
    y1_intercept_hi = outpar[1] - (err[1]);
    y2_intercept_hi = outpar[1] + (err[1]);
    
    x1_slope_lo = outpar[0] - (err[0]);
    x2_slope_lo = outpar[0] + (err[0]);
    y1_slope_lo = outpar[1] - (err[1]);
    y2_slope_lo = outpar[1] - (err[1]);
    
    x1_slope_hi = outpar[0] - (err[0]);
    x2_slope_hi = outpar[0] + (err[0]);
    y1_slope_hi = outpar[1] + (err[1]);
    y2_slope_hi = outpar[1] + (err[1]);

    
    TLine * l_intercept_lo = new TLine(x1_intercept_lo,y1_intercept_lo,x2_intercept_lo,y2_intercept_lo);
    l_intercept_lo->SetLineStyle(2);
    l_intercept_lo->SetLineColor(2);
    l_intercept_lo->Draw();
    
    TLine * l_intercept_hi = new TLine(x1_intercept_hi,y1_intercept_hi,x2_intercept_hi,y2_intercept_hi);
    l_intercept_hi->SetLineStyle(2);
    l_intercept_hi->SetLineColor(2);
    l_intercept_hi->Draw();
    
    
    TLine * l_slope_lo = new TLine(x1_slope_lo,y1_slope_lo,x2_slope_lo,y2_slope_lo);
    l_slope_lo->SetLineStyle(2);
    l_slope_lo->SetLineColor(2);
    l_slope_lo->Draw();
    
    TLine * l_slope_hi = new TLine(x1_slope_hi,y1_slope_hi,x2_slope_hi,y2_slope_hi);
    l_slope_hi->SetLineStyle(2);
    l_slope_hi->SetLineColor(2);
    l_slope_hi->Draw();
 
    c3->Write();
    c3->SaveAs("contour.png");

    TFile * fc = new TFile("chi2.root","RECREATE");
    g2->Write();
    fc->Write();
    
}



