#Function that takes dat file  and returns filled TH1F
from ROOT import TH1F
from array import *
import numpy as np

x_bins = []
y_bins_data = []
y_bins_data_unc = []
y_bins_mc = []


def make_histo_mykola_format(dat_file):
        f = open(dat_file, 'r')
        x_up = []
        bins = []
        heights = []
        stat_unc = []
        sys_unc = []
        bin_uncs =[]
        for line in f:
            if str(line.split("$")[0]) != " ":
                x = float(line.split("$")[0])
                y_data = float(line.split("$")[3])
                y_data_unc = float(line.split("$")[6])
                y_mc = float(line.split("$")[7])
                x_up.append(float(line.split("$")[2]))
                x_bins.append(x)
                y_bins_data.append(y_data)
                y_bins_data_unc.append(y_data_unc)
                y_bins_mc.append(y_mc)
        for i in range(0,len(x_up)-1):
                bin_uncs.append( (float(x_bins[i+1])  -   float(x_bins[i])) /2 )
    
        return x_bins, bin_uncs, y_bins_data, y_bins_data_unc, y_bins_mc

def make_histo_otto_format(dat_file, THEORY):
    pred_line = -1
    if THEORY == "Mitov":
        pred_line = 5
    elif THEORY == "PWHGP8":
        pred_line = 6
    f = open(dat_file, 'r')
    bin_ctrs = []
    bin_uncs = []
    heights = []
    stat_unc = []
    sys_unc = []
    line_number = 0
    for line in f:
        split_line = line.translate(None, ' "{}[];\n')
        split_line = split_line.split(",")

        if (line_number == 0):
            x_bins = line.split(",")
        if (line_number == 1):
            y_unc = line.split(",")
        if (line_number == pred_line):
            ratio = line.split(",")

        print "OTTO: = " + str(split_line)
        line_number = line_number + 1
    for i in range(0,len(x_bins)-1):
        bin_ctrs.append( float(x_bins[i])  +   ((float(x_bins[i+1])  -  float(x_bins[i]))/2 ))
        bin_uncs.append( (float(x_bins[i+1])  -   float(x_bins[i])) /2 )
    for i in y_unc:
        stat_unc.append(float(i))
    for i in ratio:
        heights.append(1.0/float(i))
    
    return bin_ctrs, bin_uncs, heights, stat_unc



def make_histo_yoda(dat_file):
    f = open(dat_file, 'r')
    x = []
    x_uncs = []
    y = []
    y_uncs = []
    line_number = 0
    for line in f:
        split_line = line.translate(None, '\n')
        split_line = split_line.split("    ")
        if ((len(split_line) > 2) & (split_line[0] != "# xval") ):
            x.append(float(split_line[0]))
            x_uncs.append(float(split_line[0])-float(split_line[1]))
            y.append(1.0 / float(split_line[3]))
            y_uncs.append(float(split_line[4])*(1.0 / float(split_line[3])))
    print str(x)
    print str(x_uncs)
    print str(y)
    print str(y_uncs)
    return x, x_uncs, y, y_uncs


























