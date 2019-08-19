import sys

if __name__ == "__main__":
    
    pretype = sys.argv[1]

    for line in sys.stdin:
        
        row = line.strip("\r\n").split("\t")
        
        if pretype == "medium":
            
            # medium format
            # rname, strand1, chrom1, pos1, frag1, strand2, chrom2, pos2, frag2, mapq1, mapq2
            
            if row[0] < row[8]:
                print " ".join(map(str, [row[3],
                                         row[4], row[0], row[2], row[7],
                                         row[12], row[8], row[10], row[15],
                                         row[6], row[14]]))
            else:
                print " ".join(map(str, [row[3],
                                         row[12], row[8], row[10], row[15],
                                         row[4], row[0], row[2], row[7],
                                         row[14], row[6]]))
