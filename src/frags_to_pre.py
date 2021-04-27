import sys

if __name__ == "__main__":
    
    pretype = sys.argv[1]
    assay = sys.argv[2]
    
    for line in sys.stdin:
        
        row = line.strip("\r\n").split("\t")
        
        if pretype == "medium":

            if assay == "BL-CAP-C":
                frag1 = 0
                frag2 = 1
            else:
                frag1 = row[7]
                frag2 = row[15]
            
            # medium format
            # rname, strand1, chrom1, pos1, frag1, strand2, chrom2, pos2, frag2, mapq1, mapq2
            
            if row[0] < row[8]:
                print " ".join(map(str, [row[3],
                                         row[4], row[0], row[2], frag1,
                                         row[12], row[8], row[10], frag2,
                                         row[6], row[14]]))
            else:
                print " ".join(map(str, [row[3],
                                         row[12], row[8], row[10], frag2,
                                         row[4], row[0], row[2], frag1,
                                         row[14], row[6]]))
