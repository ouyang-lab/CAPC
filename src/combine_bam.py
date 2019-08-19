import sys
import pysam

if __name__ == "__main__":

    """
    f_bam1 = sys.argv[1]
    f_bam2 = sys.argv[2]
    
    samy1 = pysam.AlignmentFile(f_bam1, "rb")
    samy2 = pysam.AlignmentFile(f_bam2, "rb")
    
    while True:
        try:
            read1 = samy1.next()
            read2 = samy2.next()

            print read1
            print read2
            assert False
            
            if read1.is_unmapped or read2.is_unmapped:
                continue

            if read1.query_name != read2.query_name:
                sys.stderr.write("#ERROR")
                assert False
            else:

                diff1 = read1.reference_start - read2.reference_start + 1
                diff2 = read2.reference_start - read1.reference_start + 1

                quals1 = "".join(map(str, read1.query_qualities))
                quals2 = "".join(map(str, read2.query_qualities))
                
                print "\t".join(map(str, [read1.query_name, "FLAG1",
                                          read1.reference_name, read1.reference_start + 1,
                                          read1.mapping_quality, read1.cigarstring,
                                          read2.reference_name, read2.reference_start + 1, diff1,
                                          read1.query_sequence, quals1,
                                          read1.opt]))
            
        except StopIteration:
            break

    samy1.close()
    samy2.close()
    """
    
    for no, line in enumerate(sys.stdin):
        
        if no%2 == 0:
            read1 = line.strip("\r\n").split("\t")
            continue
        elif no%2 == 1:
            read2 = line.strip("\r\n").split("\t")

        if read1[2] == "*" or read2[2] == "*":
            continue
            
        if read1[0] != read2[0]:
            
            sys.stderr.write("#ERROR\n")
            sys.stderr.write(read1[0] + "\n")
            sys.stderr.write(read2[0] + "\n")
            
            assert False
        else:

            # Paired: 1
            # Properly Paired: 2
            # Forward: 16
            # Reverse: 32
            # Read1: 64
            # Read2: 128

            if int(read1[1]) == 0:
                flag1 = 32 + 1 + 64
            elif int(read1[1]) == 16:
                flag1 = 16 + 1 + 64
            else:
                print int(read1[1])
                assert False
                
            if int(read2[1]) == 0:
                flag2 = 32 + 1 + 128
            elif int(read2[1]) == 16:
                flag2 = 16 + 1 + 128
            else:
                print int(read2[1])
                assert False
            
            if read1[2] == read2[2]:
                chrom1 = read1[2]
                chrom2 = "="
                
                diff1 = int(read2[3])-int(read1[3])
                diff2 = int(read1[3])-int(read2[3])

                # READ 1
                print "\t".join([read1[0], str(flag1),
                                 chrom1, read1[3], read1[4], read1[5],
                                 chrom2, read2[3], str(diff1),
                                 read1[9], read1[10], "\t".join(read1[11:])])

                # READ 2
                print "\t".join([read2[0], str(flag2),
                                 chrom1, read2[3], read2[4], read2[5],
                                 chrom2, read1[3], str(diff2),
                                 read2[9], read2[10], "\t".join(read2[11:])])

                
            else:
                chrom1 = read1[2]
                chrom2 = read2[2]
                diff1 = 0
                diff2 = 0
            
                # READ 1
                print "\t".join([read1[0], str(flag1),
                                 chrom1, read1[3], read1[4], read1[5],
                                 chrom2, read2[3], str(diff1),
                                 read1[9], read1[10], "\t".join(read1[11:])])

                # READ 2
                print "\t".join([read2[0], str(flag2),
                                 chrom2, read2[3], read2[4], read2[5],
                                 chrom1, read1[3], str(diff2),
                                 read2[9], read2[10], "\t".join(read2[11:])])
    
