import ROOT

def read_pred(dat_file, THEORY):
    f = open(dat_file, 'r')
    bins = []
    heights = []
    low_unc = []
    hi_unc = []
    line_number  = 0
    read = 0
    lines_read = 0
    pred = []
        
    for line in f:
        if (THEORY in line.split(" ")):
            #print "THEORY: " + THEORY
            read = 1
        if ((read == 1) & (lines_read < 9)):
            #print "READING THEORY PRED.: LINE NUMBER  " +   str(lines_read) + "  " + str(line)
            pred.append(line)
            lines_read = lines_read + 1
            """
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
                    """
    pred_line = (pred[5].split("="))[1]
    pred_line = pred_line.translate(None, ' {};\n')
    pred_line = pred_line.split(",")
    map(float, pred_line)
    for i in pred_line:
        heights.append(float(i))
    print "read_theory "
    print heights
    return heights
