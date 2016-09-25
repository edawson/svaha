import sys
import argparse

def header():
    return "\t".join(["##CHROM", "POS", "ID", "REF", "ALT", "FILTER", "INFO"])


def format_variant(line):
    ret = {}
    info = {}
    f_mat = {}
## Need a bunch of things
## sequence (chromosome/contig, should capitalize it)
## position (start position)
## ID: I guess a '.'
## Ref: N
## Alt: <TYPE>
## filter: PASS
## INFO:
##      SVTYPE
##      SVLEN
## FORMAT
#SAMPLES
   

    splits = line.strip().split("\t")
    if "inv" in splits[]:
        pass
    elif "del" in splits[]:
        if "intra" in splits[]:
            info["SVTYPE"] = "DEL"
            ret["CHROM"] = splits[0]
            ret["POS"] = splits[1]
            ret["ID"] = splits[2]
            ret["REF"] = splits[3]
            ret["ALT"] = splits[4]

            infostr = splits[
            info["SVLEN"] =

        else:
            pass
    elif "ins" in splits[]:
        pass
    return ret

def line2_vcf(d):
    ret = []
    ret.append(d["CHROM"])
    ret.append(d["POS"])
    ret.append(d["ID"])
    ret.append(d["REF"])
    ret.append(d["ALT"])

    infostr = ""
    for i in ret["INFO"]:
        infostr += i
        infostr += "="
        infostr += ret["INFO"][i]

    ret.append(infostr)

    return "\t".join(ret)

if __name__ == "__main__":
    
    print header()
    with open(sys.argv[1], "r") as fi:
        for line in fi:
            splits = line.strip().split("t")
            print line2vcf(splits)
