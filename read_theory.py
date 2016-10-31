import ROOT

def read_pred(dat_file):
    f = open(dat_file, 'r')
    bins = []
    heights = []
    low_unc = []
    hi_unc = []
    line_number  = 0
        
    for line in f:
        print "reading"
        print line.split("$")

        if (line_number == 0):
            for val in line.split("$"):
                bins.append(float(val))
        elif line_number == 1:
            for val in line.split("$"):
                heights.append(float(val))
        elif line_number == 2:
            for val in line.split("$"):
                low_unc.append(float(val))
        elif line_number == 3:
            for val in line.split("$"):
                hi_unc.append(float(val))


#  print "reading theory = = " + str(bins) + " " + str(heights) + "  " + str(low_unc) + " " + str(hi_unc)
    return  heights
