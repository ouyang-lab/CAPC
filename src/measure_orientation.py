import sys
import gzip

if __name__ == "__main__":

    # ============================================================
    # This performs strand orientation analysis on the pre files
    # ============================================================
    #
    # python measure_orientation.py <pre> <out-prefix>
    #
    
    f_pre = sys.argv[1]  # make sure it's filtered
    f_out = sys.argv[2]  # output prefix
    mapq_threshold = 1
    frag_threshold = 0
    
    f_FR = "/".join([f_out, "FR.txt.gz"])
    f_RF = "/".join([f_out, "RF.txt.gz"])
    f_FF = "/".join([f_out, "FF.txt.gz"])
    f_RR = "/".join([f_out, "RR.txt.gz"])

    FR = {}
    RF = {}
    FF = {}
    RR = {}

    total = 0
    with gzip.open(f_pre, "rb") as f:
        for line in f:
            row = line.strip("\r\n").split(" ")
            (strand1, chrom1, pos1, frag1,
             strand2, chrom2, pos2, frag2,
             mapq1, mapq2) = (int(row[1]), row[2],
                              int(row[3]), int(row[4]),
                              int(row[5]), row[6],
                              int(row[7]), int(row[8]),
                              int(row[9]), int(row[10]))
            
            # intra-chromosome, unligated, religated and MAPQ>=30
            if (chrom1 == chrom2 and abs(frag1-frag2) > frag_threshold and
                mapq1 >= mapq_threshold and mapq2 >= mapq_threshold):
                
                dist = abs(pos1-pos2)+1
                
                if strand1 == strand2:
                    if strand1 == 0 and strand2 == 0:  # FF
                        try:
                            FF[dist] += 1
                        except KeyError:
                            FF[dist] = 1
                    elif strand1 == 1 and strand2 == 1:  # RR
                        try:
                            RR[dist] += 1
                        except KeyError:
                            RR[dist] = 1
                else:
                    if strand1 == 0 and pos1 < pos2:  # FR
                        try:
                            FR[dist] += 1
                        except KeyError:
                            FR[dist] = 1
                    elif strand1 == 1 and pos1 > pos2:  # FR
                        try:
                            FR[dist] += 1
                        except KeyError:
                            FR[dist] = 1
                    else:  # RF
                        try:
                            RF[dist] += 1
                        except KeyError:
                            RF[dist] = 1

                total += 1
    
    # Output
    o_FR = gzip.open(f_FR, "wb")
    o_RF = gzip.open(f_RF, "wb")
    o_FF = gzip.open(f_FF, "wb")
    o_RR = gzip.open(f_RR, "wb")
    
    for k, v in sorted(FR.items()):
        o_FR.write("%s\t%s\t%s\n" % (k, v, float(v)/float(total)))

    for k, v in sorted(RF.items()):
        o_RF.write("%s\t%s\t%s\n" % (k, v, float(v)/float(total)))

    for k, v in sorted(FF.items()):
        o_FF.write("%s\t%s\t%s\n" % (k, v, float(v)/float(total)))

    for k, v in sorted(RR.items()):
        o_RR.write("%s\t%s\t%s\n" % (k, v, float(v)/float(total)))
    
    o_FR.close()
    o_RF.close()
    o_FF.close()
    o_RR.close()
    
