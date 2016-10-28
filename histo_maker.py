#Function that takes dat file  and returns filled TH1F
from ROOT import TH1F
from array import *
import numpy as np

x_bins = []
y_bins_data = []
y_bins_data_unc = []
y_bins_mc = []


def make_histo(dat_file):
	f = open(dat_file, 'r')
	bins = []
	heights = []
	stat_unc = []
	sys_unc = []
        for line in f:
            print "y unc (rel) = " + str(line.split("$")[6])
            if str(line.split("$")[0]) != " ":
                x = float(line.split("$")[0])
                y_data = float(line.split("$")[3])
                y_data_unc = float(line.split("$")[6])
                y_mc = float(line.split("$")[7])
                x_bins.append(x)
                y_bins_data.append(y_data)
                y_bins_data_unc.append(y_data_unc)
                y_bins_mc.append(y_mc)

        return x_bins, y_bins_data, y_bins_data_unc, y_bins_mc

"""
            if line.split("\t")[1] != "-":
			heights.append(float(line.split("\t")[1]))
			stat_unc.append(float(line.split("\t")[2]))
			sys_unc.append(float(line.split("\t")[3]))

	bins_ar = array('d',bins)
	heights_ar = array('d',heights)

       	print bins_ar

	histogram = TH1F(dat_file, dat_file, len(bins_ar)-1, bins_ar)


	for bin in range(0,len(heights)):
		total_unc = ((stat_unc[bin]**2)+(sys_unc[bin]**2))**(0.5)
		histogram.SetBinContent(bin, heights_ar[bin])
		histogram.SetBinError(bin, total_unc)


	return histogram

"""
