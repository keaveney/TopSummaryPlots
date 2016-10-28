from ROOT import *
import ROOT


def make_histo(file):

    f = TFile(file)
    c = TCanvas()
    
    data_x_n, data_y_n  = ROOT.Double(0), ROOT.Double(0)
    
    ratio_list =  []
    x_list = []
    
    histo_mc = f.Get("mc;1")
    graph_data = f.Get("data;1")

    print "b bins in MC histo = " + str(histo_mc.GetSize())
    print "b points in DATA histo = " + str(graph_data.GetN())

    for n in range(0,graph_data.GetN()):
        graph_data.GetPoint(n, data_x_n,data_y_n)
        mc = float(histo_mc.GetBinContent(n+1))
        print "N:" + str(n) + " Raw mc = " + str(mc) + "  Raw data = " + str(data_y_n)
        if (mc > 0.0):
            ratio = data_y_n/mc
            ratio_list.append(ratio)
            x_list.append(data_x_n)

    return x_list, ratio_list
