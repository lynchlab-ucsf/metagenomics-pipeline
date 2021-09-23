#!/usr/bin/env python
###################################
import gzip as gz, csv, sys

def GetFilelist(i):
    import os
    ## If user supplies an existing hist file, use it
    if os.path.exists(i) and not "/dev/" in i: return [i]
    ## Otherwise, read the files
    files = []
    hin = open(i, 'r')
    for line in hin: files.append(line.rstrip())
    hin.close()
    return files

def WriteCov(samples, genes, houtcsv):
    cover = {}
    houtcsv.writerow(["Gene"]+sorted(samples.keys()))    
    for gene in genes.keys():
        out = [gene]
        for sample in sorted(samples.keys()):
            try: c = samples[sample][gene]
            except KeyError: c = 0.0
            try: cover[sample].append(c)
            except KeyError: cover[sample] = [c]
            out.append(c)
        houtcsv.writerow(out)
    return cover

def WriteStat(cover, houtcsv):
    import numpy as np
    houtcsv.writerow(["Sample","MeanCoverage","MedianCoverage","MaxCoverage","MinCoverage"])
    for sample,l in cover.iteritems():
        median_cov = np.median(l)
        max_cov = np.max(l)
        min_cov = np.min(l)
        mean_cov = np.mean(l)
        houtcsv.writerow([sample,mean_cov,median_cov,max_cov,min_cov])

def CalcCov(files):
    samples = {}
    genes = {}
    cover = {}
    hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '\t')
    for f in files:
        if ".gz" in f: hin = gz.open(f)
        else: hin = open(f)
        sample=(f.split("/")[-1]).split(".")[0]
        samples[sample] = {}
        cover[sample] = []
        hincsv = csv.reader(hin, delimiter = '\t')
        gene = ""
        cov = 0.0
        for i,row in enumerate(hincsv):
            #if len(row)<8: break
            if gene != row[3]:
                if len(files)==1:
                    if i==0:houtcsv.writerow(["Gene",sample])
                    else: houtcsv.writerow([gene,cov])
                else: samples[sample][gene] = cov
                gene = row[3]
                cov = 0.0
            genes[gene] = ""
            depth = float(row[-4])
            depth_f = float(row[-1])
            cov += depth*depth_f
        hin.close()
        ## Handle last gene in the file
        if len(files)==1: houtcsv.writerow([gene,cov])
        else: samples[sample][gene] = cov
    if len(files)==1: sys.exit()
    cover = WriteCov(samples, genes, houtcsv)
    hout.close()
    hout = sys.stderr
    houtcsv = csv.writer(hout, delimiter = '\t')
    WriteStat(cover, houtcsv)
    hout.close()

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", type=str,
            help="Infile histogram(s) from bedtools coverage")
    args = parser.parse_args()
    if not args.infiles: sys.exit(parser.print_help())
    files = GetFilelist(args.infiles)
    CalcCov(files)

if __name__ == "__main__": main()
