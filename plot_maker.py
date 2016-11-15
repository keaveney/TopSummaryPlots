#Function that takes dat file  and returns filled TH1F
from ROOT import TH1F, TGraphErrors, TCanvas, TLegend, gStyle, TF1, TVirtualFitter, TLine, TLatex, TMathText
from array import *
import numpy as np
import read_theory
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

c = TCanvas()

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12

iPeriod = 0

gStyle.SetLegendBorderSize(0);
gStyle.SetLegendFillColor(0);
gStyle.SetLegendFont(42);
gStyle.SetLegendTextSize(0.05);


#Theory type (options are: 'Mitov',  'Pecjak',  'Lipka', 'Kidonakis', PWHGP8)
THEORY = "PWHGP8"


list_x_ll,list_x_unc_ll, list_y_data_ll, list_y_data_unc_ll, list_y_mc_ll  = histo_maker.make_histo_mykola_format("TOP_16_011.dat")


if THEORY != "PWHGP8":
    list_y_mc_ll  = read_theory.read_pred("NNLO_TOP_16_011.dat", THEORY)
    print "NNLO values = " +   str(list_y_mc_ll)

x_ll = np.asarray(list_x_ll)
y_data_ll = np.asarray(list_y_data_ll)
x_data_unc_ll = np.zeros(len(list_x_ll))
#x_data_unc_ll = np.asarray(list_x_unc_ll)
y_mc_ll = np.asarray(list_y_mc_ll)
print "data = " + str(y_data_ll)
print "mc = " + str(y_mc_ll)
ratio_ll = y_data_ll/y_mc_ll
y_data_unc_ll = (np.asarray(list_y_data_unc_ll)/100) * ratio_ll
gr_ll = TGraphErrors(len(list_x_ll), x_ll, ratio_ll, x_data_unc_ll, y_data_unc_ll)

if (THEORY == "Mitov"):
    list_x_ljets, list_x_unc_ljets, list_y_data_ljets, list_y_data_unc_ljets = histo_maker.make_histo_yoda("TOP_16_008_DATA_NNLO.yoda")
elif (THEORY == "PWHGP8"):
    list_x_ljets, list_x_unc_ljets, list_y_data_ljets, list_y_data_unc_ljets = histo_maker.make_histo_yoda("TOP_16_008_DATA_PWHGP8.yoda")

x_ljets = np.asarray(list_x_ljets)
ratio_ljets = np.asarray(list_y_data_ljets)

y_data_unc_ljets = np.asarray(list_y_data_unc_ljets)
x_data_unc_ljets = np.zeros(len(list_x_ljets))
#x_data_unc_ljets = np.asarray(list_x_unc_ljets)
y_data_ljets = np.asarray(list_y_data_ljets)
y_data_unc_ljets = np.asarray(y_data_unc_ljets)

#y_data_unc_ljets = np.subtract(y_data_unc_ljets , y_data_ljets)
y_data_unc_ljets = np.absolute(y_data_unc_ljets)

#y_data_unc_ljets = y_data_unc_ljets -  ratio_ljets
#y_data_unc_ljets = map(abs, y_data_unc_ljets)

gr_ljets = TGraphErrors(len(list_x_ljets), x_ljets, ratio_ljets , x_data_unc_ljets, y_data_unc_ljets)


xbins_comb  = np.concatenate([x_ll, x_ljets])
ratio_comb = np.concatenate([ratio_ll, ratio_ljets])
x_data_unc_comb = np.concatenate([x_data_unc_ll , x_data_unc_ljets])
y_data_unc_comb = np.concatenate([y_data_unc_ll , y_data_unc_ljets])


gr = TGraphErrors(len(xbins_comb), xbins_comb, ratio_comb, x_data_unc_comb, y_data_unc_comb)

#gr = TGraphErrors(len(list_x_ll), x_ll, ratio_ll, x_data_unc_ll, y_data_unc_ll)



# Fitting histogram (with predefined function):
fit_x = TF1("fit_x", "pol1", 0.0, 550.0)
fit_x.SetLineColor(1)
fit_x.SetLineWidth(3)
gr.Fit(fit_x)

#Create a TGraphErrors to hold the confidence intervals
grint = TGraphErrors(len(xbins_comb))
grint.SetTitle(" ")
for i in range(0,len(xbins_comb)):
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
grint.SetFillColorAlpha(16,0.65)

grint.GetHistogram().SetTitle(" ");
grint.GetHistogram().GetXaxis().SetTitle("Top P_{T} (GeV)")



if(THEORY == "PWHGP8" ):
    grint.GetHistogram().GetYaxis().SetTitle("DATA/Powheg+Pythia8")
elif(THEORY == "Kidonakis"):
    grint.GetHistogram().GetYaxis().SetTitle("DATA/aN^3LO")
elif(THEORY == "Mitov"):
    grint.GetHistogram().GetYaxis().SetTitle("DATA/NNLO")
elif(THEORY == "Lipka"):
    grint.GetHistogram().GetYaxis().SetTitle("DATA/approx. NNLO")
elif(THEORY == "Pecjak"):
    grint.GetHistogram().GetYaxis().SetTitle("DATA/approx. NLO + NNLL'")


gr_ll.SetMarkerStyle(22)
gr_ll.SetMarkerColor(4)
gr_ll.SetLineColor(4)

gr_ljets.SetMarkerStyle(23)
gr_ljets.SetMarkerColor(2)
gr_ljets.SetLineColor(2)

#grint.Draw("APE3")
#gr.Draw("PSAME")


grint.Draw("ALE3")
gr_ljets.Draw("PSAME")
gr_ll.Draw("PSAME")


grint.GetHistogram().GetXaxis().SetRangeUser(0.0, 800.0)
grint.GetHistogram().GetYaxis().SetRangeUser(0.3, 1.7)


line = TLine(0.0,1.0,660.0,1.0);
line.SetLineStyle(2)
line.Draw()

leg = TLegend(0.54,0.7,0.92,0.86)
leg.SetHeader("#bf{Parton level}")
leg.AddEntry(gr_ll,"dilepton","p")
leg.AddEntry(gr_ljets,"l+jets","p")
leg.AddEntry(grint,"linear fit","lf")

leg.Draw()

#construct string for fit results
#fit_string = "r = " + str(p2*100) + "\cdot( \frac{P_{T}}{100\; GeV} ) + " + str(p1)
#fit_string = "r = " + str(p2*100) +  " \\cdot {P_{T}} \\over {100\;GeV} + " + str(p1)

fit_string = "r = " + str(p2*100) + "\cdot({P_{T}\\over100\;GeV}) + "+ str(p1)

#display fit results on plot
#latex = TLatex()
#latex.SetTextSize(0.075)
#latex.SetTextAlign(13)
#latex.DrawLatex(50.0,0.6, fit_string)

l = TMathText()
l.SetTextAlign(23)
l.SetTextSize(0.04)
l.DrawMathText(250.0, 0.6, fit_string);

#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c, iPeriod, iPos)


if(THEORY == "Kidonakis"):
    c.SaveAs("Ratio_DATA_Kidonakis.png")
elif(THEORY == "Lipka"):
    c.SaveAs("Ratio_DATA_Lipka.png")
elif(THEORY == "Pecjak"):
    c.SaveAs("Ratio_DATA_Pecjak.png")
elif(THEORY == "Mitov"):
    c.SaveAs("Ratio_DATA_Mitov.png")
elif(THEORY == "PWHGP8"):
    c.SaveAs("Ratio_DATA_PWHGP8.png")








