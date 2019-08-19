import sys
import pysam

if __name__ == "__main__":

    f_bam = sys.argv[1]
    f_read1 = sys.argv[2]
    f_read2 = sys.argv[3]

    o_read1 = open(f_read1, "w")
    o_read2 = open(f_read2, "w")
    
    current_qname = None
    current_count = 0
    current_toskip = False
    read1 = []
    read2 = []

    counter_gt_2 = 0
    
    
    samy = pysam.AlignmentFile(f_bam, "rb")
    for read in samy:
        
        qname = read.query_name
        if read.is_unmapped:
            chrom = "*"
        else:
            chrom = read.reference_name
        
        start = read.reference_start  # 0-based
        end = read.reference_end  # 1-based
        strand = 1 if read.is_reverse else 0
        rp = 0 if read.is_read1 else 1
        mapq = read.mapping_quality
        
        if current_qname != qname:

            if current_qname != None:
                if not current_toskip:
                    if current_count == 2:
                        o_read1.write("\t".join(map(str, read1[0])) + "\n")
                        o_read2.write("\t".join(map(str, read2[0])) + "\n")
                    elif current_count > 2:
                        counter_gt_2 += 1
            
            current_qname = qname
            current_count = 0
            current_toskip = False
            read1 = []
            read2 = []
        
        if read.is_read1:

            if read.is_unmapped:
                current_toskip = True
            else:
                if not read.is_reverse:
                    read1.append([chrom, start, start+1, qname, strand, rp, mapq])
                else:
                    read1.append([chrom, end-1, end, qname, strand, rp, mapq])
                    
        elif read.is_read2:

            if read.is_unmapped:
                current_toskip = True
            else:
                if not read.is_reverse:
                    read2.append([chrom, start, start+1, qname, strand, rp, mapq])
                else:
                    read2.append([chrom, end-1, end, qname, strand, rp, mapq])
        
        current_count += 1

    sys.stderr.write("Greater than 2: %s\n" % counter_gt_2)

    o_read1.close()
    o_read2.close()
