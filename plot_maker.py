#Function that takes dat file  and returns filled TH1F
from ROOT import TH1F, TGraphErrors, TCanvas, TLegend, gStyle, TF1, TVirtualFitter, TLine, TLatex, TMathText
from array import *
import numpy as np
import histo_maker
import histo_maker_fromROOT
import tdrstyle
import CMS_lumi


#set the tdr style
tdrstyle.setTDRStyle()


#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_13TeV = "2.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "(13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)


fromROOT = 0
fromDAT  = 1

c = TCanvas()

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12

iPeriod = 0


gStyle.SetLegendBorderSize(0);
gStyle.SetLegendFillColor(0);
gStyle.SetLegendFont(42);
gStyle.SetLegendTextSize(0.05);

if fromDAT == 1:
    print "making ratio from DAT file"
    list_x, list_y_data, list_y_data_unc, list_y_mc  = histo_maker.make_histo("TOP_16_011.dat")
    x = np.asarray(list_x)
    y_data = np.asarray(list_y_data)
    x_data_unc = np.zeros(len(list_x))
    y_mc = np.asarray(list_y_mc)
    ratio = y_data/y_mc
    y_data_unc = (np.asarray(list_y_data_unc)/100) * ratio
    print "data = " + str(y_data)
    print "mc = " + str(y_mc)
    gr = TGraphErrors(len(list_x), x, ratio, x_data_unc, y_data_unc)
    

    # Fitting histogram (with predefined function):
    fit_x = TF1("fit_x", "pol1", 0.0, 550.0)
    fit_x.SetLineColor(1)
    fit_x.SetLineWidth(3)
    gr.Fit("fit_x")

    #Create a TGraphErrors to hold the confidence intervals
    grint = TGraphErrors(len(list_x))
    grint.SetTitle(" ")
    for i in range(0,len(list_x)):
       grint.SetPoint(i, gr.GetX()[i], 0)
       #Compute the confidence intervals at the x points of the created graph
       #(TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint);
       #Now the "grint" graph contains function values as its y-coordinates
       #and confidence intervals as the errors on these coordinates
       #Draw the graph, the function and the confidence intervals
       
       
    #fit = hist->GetFunction(function_name);
    chi2 = fit_x.GetChisquare()
    # value of the first parameter
    p1 = fit_x.GetParameter(0)
    p1 = round(p1,4)
    # error of the first parameter
    e1 = fit_x.GetParError(0)
    e1 = round(e1,4)

    # value of the first parameter
    p2 = fit_x.GetParameter(1)
    p2 = round(p2,4)

    # error of the first parameter
    e2 = fit_x.GetParError(1)
    e2 = round(e2,4)


    print "Chi2 = " + str(chi2)
    print "P1 = " + str(p1)
    print "E1 = " + str(e1)


    TVirtualFitter.GetFitter().GetConfidenceIntervals(grint)
    grint.SetLineColor(1)
    grint.SetFillColor(16)

    gr.GetHistogram().SetTitle(" ");
    grint.GetHistogram().GetXaxis().SetTitle("Top P_{T} (GeV)")
    grint.GetHistogram().GetYaxis().SetTitle("DATA/Powheg+Pythia8")
    grint.GetHistogram().GetYaxis().SetRangeUser(0.4, 1.6)

    gr.SetMarkerStyle(22)
    gr.SetMarkerColor(4)
    gr.SetLineColor(4)

    grint.Draw("APE3")
    gr.Draw("PSAME")

    line = TLine(0.0,1.0,520.0,1.0);
    line.SetLineStyle(2)
    line.Draw()

    leg = TLegend(0.48,0.6,0.88,0.8)
    #leg.SetHeader("Parton level")
    leg.AddEntry(gr,"dilepton","p")
    leg.Draw()

    #construct string for fit results
    #    fit_string = "r = " + str(p2*100) + "\cdot( \frac{P_{T}}{100\; GeV} ) + " + str(p1)
    # fit_string = "r = " + str(p2*100) +  " \\cdot {P_{T}} \\over {100\;GeV} + " + str(p1)

    fit_string = "r = " + str(p2*100) + "\cdot({P_{T}\\over100\;GeV}) + "+ str(p1)

    #display fit results on plot
    # latex = TLatex()
    # latex.SetTextSize(0.075)
    # latex.SetTextAlign(13)
    # latex.DrawLatex(50.0,0.6, fit_string)

    l = TMathText()
    l.SetTextAlign(23)
    l.SetTextSize(0.04)
    l.DrawMathText(170.0, 0.6, fit_string);



    #draw the lumi text on the canvas
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)


if fromROOT == 1:

    list_x, list_y, norm = histo_maker_fromROOT.make_histo("TOP_16_011.root")
    print list_x
    print list_y

    x = np.asarray(list_x)
    y = np.asarray(list_y)

    gr = TGraphErrors(len(list_x), x, y)

    gr.Draw("AP")

c.SaveAs("Ratio.pdf")
