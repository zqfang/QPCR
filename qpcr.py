#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import math
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python Qpcr.py --input date.txt
input: 
Well	Fluor	Target	Content	Sample	Cq
B02	SYBR	ubq	Unkn-01	nb	20.19
B03	SYBR	ubq	Unkn-01	nb	20.24
B04	SYBR	ubq	Unkn-01	nb	20.30
B05	SYBR	ubq	Unkn-01	nb	20.29
B06	SYBR	ubq	Unkn-02	heg4	18.75
B08	SYBR	ubq	Unkn-02	heg4	18.67
B09	SYBR	ubq	Unkn-02	heg4	18.85
C02	SYBR	ef1	Unkn-03	nb	18.14
C03	SYBR	ef1	Unkn-03	nb	18.23
C05	SYBR	ef1	Unkn-03	nb	18.16
C06	SYBR	ef1	Unkn-04	heg4	16.62
C07	SYBR	ef1	Unkn-04	heg4	16.60
C09	SYBR	ef1	Unkn-04	heg4	16.77
    '''
    print message

def sqr(x):
    return x*x

def readinput(infile, outfile):
    ctl_gene = 'ef1'
    ctl_sample = 'nb'
    '''Raw data input'''
    data = defaultdict(lambda : defaultdict(lambda: list()))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('Well'): 
                unit = re.split(r'\t',line)
                data[unit[2]][unit[4]].append(float(unit[5]))
    '''Mean and std'''
    data_sum = defaultdict(lambda : defaultdict(lambda: list())) 
    for gene in data.keys():
        for sample in data[gene].keys():
            #print gene, sample, np.mean(data[gene][sample]), np.std(data[gene][sample]), data[gene][sample]
            data_sum[gene][sample] = [np.mean(data[gene][sample]), np.std(data[gene][sample])]
            #print gene, sample, data_sum[gene][sample]
            

    #DeltaCt
    data_expr = defaultdict(lambda : defaultdict(lambda: list()))
    for gene in sorted(data_sum.keys()):
        if gene not in ctl_gene:
            for sample in data_sum[gene].keys():
                DeltaCt = float(data_sum[gene][sample][0]) - float(data_sum[ctl_gene][sample][0])
                S_std   = math.sqrt(sqr(float(data_sum[gene][sample][1])) + sqr(float(data_sum[ctl_gene][sample][1])))
                data_expr[gene][sample] = [DeltaCt, S_std]

    #2^(-DeltaDeltaCt): DeltaDeltaCt2
    samples = defaultdict(int)
    data_final = defaultdict(lambda : defaultdict(lambda: list()))
    for gene in sorted(data_expr.keys()):
        if gene not in ctl_gene:
            for sample in data_expr[gene].keys():
                #print gene, sample
                samples[sample] = 1
                if sample in ctl_sample:
                    DeltaDeltaCt = 0.00
                    DeltaDeltaCt2 = 1.00
                    Range_botton = 2**(-DeltaDeltaCt-data_expr[gene][sample][1])
                    Range_top    = 2**(-DeltaDeltaCt+data_expr[gene][sample][1])
                    STDEV        = ((DeltaDeltaCt2-Range_botton)+(Range_top-DeltaDeltaCt2))/2
                    data_final[gene][sample] = [DeltaDeltaCt2, STDEV]
                else:
                    DeltaDeltaCt = data_expr[gene][sample][0] - data_expr[gene][ctl_sample][0]
                    DeltaDeltaCt2 = 2**(-DeltaDeltaCt)
                    Range_botton = 2**(-DeltaDeltaCt-data_expr[gene][sample][1])
                    Range_top    = 2**(-DeltaDeltaCt+data_expr[gene][sample][1])
                    STDEV        = ((DeltaDeltaCt2-Range_botton)+(Range_top-DeltaDeltaCt2))/2 
                    data_final[gene][sample] = [DeltaDeltaCt2, STDEV]
    headers = ['Gene']
    for s in sorted(samples.keys()):
        headers.append(s)
        headers.append('STDEV')
    
    ofile = open(outfile, 'w')
    header = '\t'.join(headers)   
    #header = '\t'.join(samples)
    print >> ofile, header
    for gene in sorted(data_final.keys()):
        expr = []
        for sample in sorted(data_final[gene].keys()):
            expr.append(data_final[gene][sample][0])
            expr.append(data_final[gene][sample][1])
        print >> ofile, '%s\t%s\t%s\t%s\t%s' % (gene, expr[0], expr[1], expr[2], expr[3])
    ofile.close()
    R_cmd(outfile)

def R_cmd(exprsumfile):
    R = '''
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
pdf("Qpcr.pdf")
data =read.table("''' + exprsumfile  + '''", header = T)
expr = rbind(data[,2], data[,4])
std = rbind(data[,3], data[,5])
barx <- barplot(expr, beside=TRUE, col=c("blue", "orange"), ylim=c(0,2), border=F, axis.lty=1, xlab="Gene", ylab="Relative Expression Level")
error.bar(barx, expr, std)
axis(1,c(0.9,max(barx)+0.6),line=0,labels=c("",""))
text(barx[1,]+0.5,rep(-0.07,6),offset=2,labels=data$Gene,srt=0,xpd=TRUE)
legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))
dev.off()
''' 
    infile = open ('Qpcr.R', 'w')
    print >> infile, R
    infile.close() 
    cmd = 'cat Qpcr.R | R --slave'
    os.system(cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if args.output is None:
        args.output = 'Qpcr.sum'
    readinput(args.input, args.output)

if __name__ == '__main__':
    main()
