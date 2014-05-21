__author__ = 'mikeknowles'
import argparse
from Bio import SeqIO
import time
import re
import os

def seperator(i,path):
    os.system("rm " + path + "/*")
    start = time.strftime("%H:%M:%S")
    counter = {}
    print "[%s] File location is: %s" % (start, i)
    for seq in SeqIO.parse(i,"fasta"):
        current = time.strftime("%H:%M:%S")
        # print seq.description
        try:
            gene = re.search('\((\w+)\)$', seq.description)
            genename = gene.group(1)
        except:
            gene = re.search('\((\w+)\)\s\d+-\d+$', seq.description)
            genename = gene.group(1)
        genepath = "%s/%s.fasta" % (os.path.expanduser(path), genename)
        genefile = open(genepath, 'a')
        SeqIO.write(seq, genefile, "fasta")
        genefile.close()
        try:
            counter[genename] = counter.get(genename, 0) + 1
        except:
            counter[genename] = 1
        print "[%s] %s sequence(s) printed to %s.fasta" % (current, counter[genename], genename)
    print path

parser = argparse.ArgumentParser(description='Convert fasta file to genes fasta')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
parser.add_argument('-i', '--input', required=True, help='Specify input fasta')
parser.add_argument('-o', '--output', required=True, help='Specify output folder')
args = vars(parser.parse_args())
seperator(args['input'], args['output'])
# seperator("/nas/Pipeline_development/AntimicrobialResistance/ARMI/db/abrg.fa","/nas/Pipeline_development/AntimicrobialResistance/ARMI/abrgdb/")

