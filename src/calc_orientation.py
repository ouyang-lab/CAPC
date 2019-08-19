import sys
import gzip
import numpy as np

if __name__ == "__main__":

    f_names = sys.argv[1:]
    max_value = 100000
    bin_size = 50
    threshold = 0.01
    
    data = []

    total_bins = (max_value/bin_size)+1
    
    for no, f_name in enumerate(f_names):

        #prefix = f_name.split("/")[-1].replace(".txt.gz", "")

        d = np.zeros(total_bins)
        
        with gzip.open(f_name, "rb") as f:
            for line in f:
                row = line.strip("\r\n").split("\t")
                size, count = (int(row[0]), int(row[1]))

                if size < max_value:
                    s = size/bin_size
                    d[s] += count
                else:
                    d[max_value/bin_size] += count

        d = d[::-1].cumsum()
        
        data.append(d)

    data = np.array(data)

    current_size = max_value
    for no, d in enumerate(data.T):
        p = d/d.sum()
        if np.all(abs(p-0.25)<=threshold):
            current_size = (total_bins-no)*bin_size
        else:
            break

    print "Orientation Size (+/-%s): %s" % (threshold, current_size)

    for no, d in enumerate(data.T):
        p = d/d.sum()
        print "\t".join(map(str, [(total_bins-no)*bin_size]+p.tolist()))
    
