from __future__ import print_function
import sys
import argparse
import gzip

def parse_cutadapt(i_P1, i_P2, o_P1, o_P2, o_F1, o_F2):
    no = 0
    counter = 0
    counter_actual_length = 0
    counter_min_length = 0
    counter_same_merge = 0
    while True:
        try:
            l1 = i_P1.next()
            l2 = i_P2.next()

            m = no % 4
            if m == 0:
                i1 = l1.rstrip("\r\n").split()[0]
                i2 = l2.rstrip("\r\n").split()[0]
                if i1 != i2:
                    print(i1)
                    print(i2)
                    assert False
            elif m == 1:
                s1 = l1.rstrip("\r\n")
                s2 = l2.rstrip("\r\n")
            elif m == 3:
                q1 = l1.rstrip("\r\n")
                q2 = l2.rstrip("\r\n")

                counter += 1
                
                if ((len(s1) != p_actual_length) or
                    (len(s2) != p_actual_length)):  # test it
                    
                    counter_actual_length += 1
                    
                    if ((len(s1) >= p_min_length) and
                        (len(s2) >= p_min_length)):

                        counter_min_length += 1

                        if s1 != s2:

                            counter_same_merge += 1
                            
                            print("\n".join(map(str, [i1, s1, "+", q1])), file=o_P1)
                            print("\n".join(map(str, [i2, s2, "+", q2])), file=o_P2)

                        else:
                            if "MNOLIG" in p_allow:
                                print("\n".join(map(str, [i1, s1, "+", q1])), file=o_F1)
                                print("\n".join(map(str, [i2, s2, "+", q2])), file=o_F2)

                else:
                    if "MNOLIG" in p_allow:
                        print("\n".join(map(str, [i1, s1, "+", q1])), file=o_F1)
                        print("\n".join(map(str, [i2, s2, "+", q2])), file=o_F2)
                    
            no += 1
        except StopIteration:
            break

    return (counter, counter_actual_length, counter_min_length, counter_same_merge)

def parse_masked_cutadapt(o_OK1, o_OK2, o_S1, o_S2, o_N1, o_N2):

    counter_total = 0
    counter_kept = 0
    counter_no_ligation = 0
    counter_multi_split = 0
    counter_too_short = 0
    counter_OK = 0
    counter_unsure = 0
    
    for no, line in enumerate(sys.stdin):
        
        row = line.strip("\r\n").split("\t")
        
        if (no % 4) == 0:
            R1 = row[0]
            R2 = row[2]
        elif (no % 4) == 3:
            if use and curr_f_handler:
                q_G1 = row[0][l_G1:lr-r_G1]
                q_A1 = row[1][l_A1:lr-r_A1]
                q_G2 = row[2][l_G2:lr-r_G2]
                q_A2 = row[3][l_A2:lr-r_A2]
                Q = [q_G1, q_A1, q_G2, q_A2]
                
                print("%s" % Q[use[0]], file=curr_f_handler[0])
                print("%s" % Q[use[1]], file=curr_f_handler[1])
                
        elif (no % 4) == 1:

            counter_total += 1
            
            use = None
            curr_f_handler = None
            lr = len(row[0])
            
            (l_G1, l_A1, l_G2, l_A2) = (lr-len(row[0].lstrip("N")), lr-len(row[1].lstrip("N")),
                                        lr-len(row[2].lstrip("N")), lr-len(row[3].lstrip("N")))
            (r_G1, r_A1, r_G2, r_A2) = (lr-len(row[0].rstrip("N")), lr-len(row[1].rstrip("N")),
                                        lr-len(row[2].rstrip("N")), lr-len(row[3].rstrip("N")))

            G1 = row[0].strip("N")  # 5' read 1
            A1 = row[1].strip("N")  # 3' read 1
            G2 = row[2].strip("N")  # 5' read 2
            A2 = row[3].strip("N")  # 3' read 2
            
            no_of_mask = sum([ 1 if mask in seq else 0 for seq in row ])
            if no_of_mask > 0:

                d_G1, d_A1, d_G2, d_A2 = len(G1), len(A1), len(G2), len(A2)

                if (d_G1 <= p_thresh or d_A1 <= p_thresh) and (d_G2 <= p_thresh or d_A2 <= p_thresh):
                    #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                    # use difference to decide

                    counter_multi_split += 1

                elif (d_G1 <= p_thresh or d_A1 <= p_thresh) and (d_G2 > p_thresh and d_A2 > p_thresh):

                    if (d_G1 <= p_thresh and d_G1 >= p_min_length) and (d_A1 <= p_thresh and d_A1 >= p_min_length):

                        if "OK" in p_allow:
                            use = (0, 1)
                            curr_f_handler = (o_OK1, o_OK2)
                            #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                            #DEBUG print G1, A1
                            print("%s\n%s\n+" % (R1, G1), file=o_OK1)
                            print("%s\n%s\n+" % (R2, A1), file=o_OK2)

                            counter_kept += 1
                            
                        counter_OK += 1
                    else:

                        if "SHORT" in p_allow:
                            
                            curr_f_handler = (o_S1, o_S2)
                            
                            if d_G1 <= p_thresh and d_G1 > d_A1:
                                use = (0, 2)
                                #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                                #DEBUG print G1, G2
                                print("%s\n%s\n+" % (R1, G1), file=o_S1)
                                print("%s\n%s\n+" % (R2, G2), file=o_S2)
                            elif d_A1 <= p_thresh and d_G1 < d_A1:
                                use = (1, 3)
                                #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                                #DEBUG print A1, A2
                                print("%s\n%s\n+" % (R1, A1), file=o_S1)
                                print("%s\n%s\n+" % (R2, A2), file=o_S2)
                            elif d_G1 > p_thresh and d_A1 < d_G1:
                                use = (1, 3)
                                #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                                #DEBUG print A1, A2
                                print("%s\n%s\n+" % (R1, A1), file=o_S1)
                                print("%s\n%s\n+" % (R2, A2), file=o_S2)
                            elif d_A1 > p_thresh and d_G1 < d_A1:
                                use = (0, 2)
                                #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                                #DEBUG print G1, G2
                                print("%s\n%s\n+" % (R1, G1), file=o_S1)
                                print("%s\n%s\n+" % (R2, G2), file=o_S2)

                            counter_kept += 1
                            
                        counter_too_short += 1

                elif (d_G1 > p_thresh and d_A1 > p_thresh) and (d_G2 <= p_thresh or d_A2 <= p_thresh):

                    if (d_G2 <= p_thresh and d_G2 >= p_min_length) and (d_A2 <= p_thresh and d_A2 >= p_min_length):
                        
                        if "OK" in p_allow:
                            use = (2, 3)
                            curr_f_handler = (o_OK1, o_OK2)
                            #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                            #DEBUG G2, A2
                            print("%s\n%s\n+" % (R1, G2), file=o_OK1)
                            print("%s\n%s\n+" % (R2, A2), file=o_OK2)

                            counter_kept += 1
                            
                        counter_OK += 1
                    else:

                        if "SHORT" in p_allow:

                            curr_f_handler = (o_S1, o_S2)
                            
                            if d_G2 <= p_thresh and d_G2 > d_A2:
                                use = (0, 2)
                                #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                                #DEBUG print G1, G2
                                print("%s\n%s\n+" % (R1, G1), file=o_S1)
                                print("%s\n%s\n+" % (R2, G2), file=o_S2)
                            elif d_A2 <= p_thresh and d_G2 < d_A2:
                                use = (1, 3)
                                #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                                #DEBUG print A1, A2
                                print("%s\n%s\n+" % (R1, A1), file=o_S1)
                                print("%s\n%s\n+" % (R2, A2), file=o_S2)
                            elif d_G2 > p_thresh and d_A2 < d_G2:
                                use = (1, 3)
                                #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                                #DEBUG print A1, A2
                                print("%s\n%s\n+" % (R1, A1), file=o_S1)
                                print("%s\n%s\n+" % (R2, A2), file=o_S2)
                            elif d_A2 > p_thresh and d_G2 < d_A2:
                                use = (0, 2)
                                #DEBUG print "\t".join(map(str, [len(G1), len(A1), len(G2), len(A2)]))
                                #DEBUG print G1, G2
                                print("%s\n%s\n+" % (R1, G1), file=o_S1)
                                print("%s\n%s\n+" % (R2, G2), file=o_S2)

                            counter_kept += 1
                            
                        counter_too_short += 1

                else:
                    counter_unsure += 1
            else:

                if "NOLIG" in p_allow:
                    use = (0, 2)
                    curr_f_handler = (o_N1, o_N2)
                    print("%s\n%s\n+" % (R1, G1), file=o_N1)
                    print("%s\n%s\n+" % (R2, G2), file=o_N2)
                    
                    counter_kept += 1
                    
                counter_no_ligation += 1
                
    return (counter_total, counter_kept, counter_no_ligation,
            counter_multi_split, counter_too_short,
            counter_OK, counter_unsure)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-p", help="app")
    parser.add_argument("-s", help="actual length", type=int)
    parser.add_argument("-m", help="minimum length", type=int, default=20)
    parser.add_argument("-n", help="match length", type=int, default=10)
    parser.add_argument("-m1", help="merge P1")
    parser.add_argument("-m2", help="merge P2")
    parser.add_argument("-o1", help="output P1")
    parser.add_argument("-o2", help="output P2")
    parser.add_argument("-prefix", default=None, help="prefix")
    parser.add_argument("-e", help="error/stats file")
    parser.add_argument("-allow", default="OK,SHORT,NOLIG", help="")
    
    args = parser.parse_args()

        
    mask = "NNNNNNNNNN"
    p_allow = args.allow.split(",")  # ["NOLIG", "OK", "SHORT"]
    p_actual_length = args.s
    p_match_length = args.n
    p_min_length = args.m
    p_thresh = p_actual_length - p_match_length
    f_prefix = args.prefix
    
    o_m_T1 = gzip.open(args.m1, "rb")
    o_m_T2 = gzip.open(args.m2, "rb")
    o_e = open(args.e, "w")
    
    if f_prefix:
        p_allow = ["OK", "SHORT", "NOLIG", "MNOLIG"]
        o_merge_nolig1 = gzip.open(f_prefix + "_MN_R1.fastq.gz", "wb")
        o_merge_nolig2 = gzip.open(f_prefix + "_MN_R2.fastq.gz", "wb")
        o_unmerge_nolig1 = gzip.open(f_prefix + "_UN_R1.fastq.gz", "wb")
        o_unmerge_nolig2 = gzip.open(f_prefix + "_UN_R2.fastq.gz", "wb")
        o_unmerge_SHORT1 = gzip.open(f_prefix + "_US_R1.fastq.gz", "wb")
        o_unmerge_SHORT2 = gzip.open(f_prefix + "_US_R2.fastq.gz", "wb")
        o_OK1 = gzip.open(f_prefix + "_OK_R1.fastq.gz", "wb")  # both merge and unmerge
        o_OK2 = gzip.open(f_prefix + "_OK_R2.fastq.gz", "wb")  # both merge and unmerge


        # Merge
        counter_merge = parse_cutadapt(o_m_T1, o_m_T2,
                                       o_OK1, o_OK2,
                                       o_merge_nolig1,
                                       o_merge_nolig2)
        # Unmerge
        (counter_total, counter_kept,
         counter_no_ligation, counter_multi_split,
         counter_too_short, counter_OK, counter_unsure) = parse_masked_cutadapt(o_OK1, o_OK2,
                                                                                o_unmerge_SHORT1,
                                                                                o_unmerge_SHORT2,
                                                                                o_unmerge_nolig1,
                                                                                o_unmerge_nolig2)
        
    else:
        o_P1 = gzip.open(args.o1, "wb")
        o_P2 = gzip.open(args.o2, "wb")

        # Merge
        counter_merge = parse_cutadapt(o_P1, o_P2,
                                       o_P1, o_P2,
                                       o_P1, o_P2)
        
        # Unmerge
        (counter_total, counter_kept,
         counter_no_ligation, counter_multi_split,
         counter_too_short, counter_OK, counter_unsure) = parse_masked_cutadapt(o_P1, o_P2,
                                                                                o_P1, o_P2,
                                                                                o_P1, o_P2)
    
    

    
    print("\t".join(["TOTAL", "KEPT", "NOLIG", "OK", "MultiSplit", "SHORT", "UNSURE"]), file=o_e)
    print("\t".join(map(str, [counter_merge[0], counter_merge[3],
                              (counter_merge[0]-counter_merge[1])+(counter_merge[2]-counter_merge[3]),
                              counter_merge[3], 0,
                              counter_merge[1]-counter_merge[2], 0])), file=o_e)
    print("\t".join(map(str, [counter_total, counter_kept,
                              counter_no_ligation, counter_OK, counter_multi_split,
                              counter_too_short, counter_unsure])), file=o_e)



    if f_prefix:
        o_merge_nolig1.close()
        o_merge_nolig2.close()
        o_unmerge_nolig1.close()
        o_unmerge_nolig2.close()
        o_unmerge_SHORT1.close()
        o_unmerge_SHORT2.close()
        o_OK1.close()
        o_OK2.close()
    else:
        o_P1.close()
        o_P2.close()

    o_m_T1.close()
    o_m_T2.close()
    o_e.close()
