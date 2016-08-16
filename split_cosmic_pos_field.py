import sys
import re

if __name__ == "__main__":

    with open(sys.argv[1], "r") as ifi:
        chrom_to_start_end = {}
        count = 1
        for line in ifi:
            line = line.strip()
            splits = line.split(":")
            if len(splits) < 2:
                continue
            CHROM = splits[0]

            splits[1] = re.sub("g.", "", splits[1])
            splits[1] = re.sub("del", "", splits[1])
            pos_l = splits[1].split("_")
            sorted_pos_l = sorted([ int( i ) for i in pos_l])

            svlen = sorted_pos_l[1] - sorted_pos_l[0]
            svlen = str(svlen)
            info = "SVLEN=" + svlen + ";SVTYPE=DEL"

            print "\t".join([CHROM, str(sorted_pos_l[0]) , str(count),
                            "N", "<DEL>", ".", ".", info, "", ""
                            ])
            count+=1

