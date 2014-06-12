__author__ = 'mikeknowles'
""" Includes threading found in examples:
http://www.troyfawkes.com/learn-python-multithreading-queues-basics/
http://www.ibm.com/developerworks/aix/library/au-threadingpython/
https://docs.python.org/2/library/threading.html
"""
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML
from threading import Thread
from Queue import Queue
from collections import defaultdict
import subprocess, os, glob, time, sys, shlex, argparse, re, csv

count = 0


def dotter():
    global count
    if count <= 80:
        sys.stdout.write('.')
        count += 1
    else:
        sys.stdout.write('\n[%s].' % (time.strftime("%H:%M:%S")))
        count = 0


def makeblastdb(dqueue):
    while True:
        #grabs fastapath from dqueue
        fastapath = dqueue.get()
        #makeblastdb if not exist
        nhr = "%s.nhr" % (fastapath)
        if not os.path.isfile(str(nhr)):
            print "makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, fastapath)
            subprocess.Popen(shlex.split("makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, fastapath)))
            dotter()
        # signals to dqueue job is done
        dqueue.task_done()

dqueue = Queue()
blastqueue = Queue()
parsequeue = Queue()
plusqueue = Queue()

def makedbthreads(fastas):
    ''' Setup and create threads for class'''
    for i in range(len(fastas)):
        threads = Thread(target=makeblastdb, args=(dqueue,))
        threads.setDaemon(True)
        threads.start()
    for fasta in fastas:
        dqueue.put(fasta)
    #wait on the dqueue until everything has been processed
    dqueue.join()

def runblast(blastqueue):
    while True:
        genome, db, out = blastqueue.get()
        blastn = NcbiblastnCommandline(query=genome, db=db, evalue=1e-40, out=out, outfmt=5)
        stdout, stderr = blastn()
        blastqueue.task_done()


def blastnthreads(fastas, genomes):
    global antidict
    for i in range(len(fastas)):
        threads = Thread(target=runblast, args=(blastqueue,))
        threads.setDaemon(True)
        threads.start()
    blastpath = {}
    for genome in genomes:
        for fasta in fastas:
            gene = re.search('\/(\w+)\.fasta', fasta)
            path = re.search('(.+)\/(.+?)\.fasta', genome)
            out = "%s/../tmp/%s.%s.xml" % (path.group(1), path.group(2), gene.group(1))
            blastpath[out] = {path.group(2): {gene.group(1)}}
            if not os.path.isfile(out):
                blastqueue.put((genome, fasta, out))
                dotter()
        #wait on the queue until everything has been processed
        blastqueue.join()
    return blastpath

plusdict = {}


def blastparser(hsp, genomes):
    global plusdict
    if hsp.identities == hsp.align_length:
        plus = '+'
    elif hsp.align_length > hsp.identities >= (hsp.align_length * 0.9):
        perid = int((float(hsp.identities) / hsp.align_length) * 100)
        plus = '(%s%%)' % (perid)
    else:
        plus = '-'
    for genome in genomes:
        for gene in genomes[genome]:
            try:
                if plus == '+':
                    plusdict[genome][gene] = {plus}
                elif plusdict[genome][gene] != {'+'}:
                    plusdict[genome][gene] = {plus}
                # elif plus == '(%s%%)' % (perid):
                #     print "%s %s %s" % (gene , genome, plus)
            except:
                if genome not in plusdict:
                    plusdict[genome] = {gene: {plus}}
                elif gene not in plusdict[genome]:
                    plusdict[genome][gene] = {plus}
                elif plus == '+':
                    plusdict[genome][gene] = {plus}

def parsethreader(xml, genomes):
    global plusdict
    dotter()
    numhsp = 0
    with open(xml) as file:
        numhsp = sum(line.count('<Hsp>') for line in file)
    if numhsp != 0:
        handle = open(xml)
        records = NCBIXML.parse(handle)
        for record in records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    blastparser(hsp, genomes)
        # parsequeue.join()
    else:
        for genome in genomes:
            for gene in genomes[genome]:
                try:
                    if gene not in plusdict[genome]:
                        plusdict[genome][gene] = {'-'}
                except:
                    plusdict[genome] = {gene: {'-'}}

def defdict(file):
    ndict = defaultdict(list)
    with open(args['tab']+file) as tab:
        for line in tab:
            tline = line.strip().split(None, 1)
            ndict[tline[0]].append(tline[1])
    return ndict

def blaster(path, targets, out):
    #modify global dotter() counter
    global count
    #retrieve markers from input
    fastas = glob.glob(path + "*.fasta")
    #retrieve genomes from input
    genomes = glob.glob(targets + "*.fasta")
    sys.stdout.write("[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S")))
    #push markers to threads
    makedbthreads(fastas)
    print "\n[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S"))
    print "[%s] Now performing BLAST database searches" % (time.strftime("%H:%M:%S"))
    blastpath = blastnthreads(fastas, genomes)
    print "[%s] Now parsing BLAST database searches" % (time.strftime("%H:%M:%S"))
    count = 80
    require = defdict("require.tab")
    antidict = defdict("typeResis.tab")
    genedict = defaultdict(dict)
    for xml in blastpath:
        parsethreader(xml, blastpath[xml])
    csvheader = 'Strain'
    row = ""
    antirow = ""
    rowcount = 0
    antihead = "Strain,Gene,\"Resistance Antibiotic\",\"Unique Resistance\"\n"
    for genomerow in plusdict:
        genedict[genomerow] = {'anti': [], 'gene': [], 'unique': []}
        row += "\n" + genomerow
        rowcount += 1
        for generow in sorted(plusdict[genomerow]):
            if rowcount <= 1:
                csvheader += ',' + generow
            tempcount = 0
            if generow.lower() in require:
                for othergene in require[generow.lower()]:
                    for propercasegene in plusdict[genomerow]:
                        if othergene == propercasegene.lower:
                            if plusdict[genomerow][propercasegene] != {'-'}:
                                tempcount += 1
            else:
                tempcount = 1
            # print tempcount, len(require[generow.lower()])
            if tempcount == len(require[generow.lower()]):
                genedict[genomerow]['gene'].append(generow)
                for antibio in antidict[generow.lower()]:
                    if antibio not in genedict[genomerow]['anti']:
                        genedict[genomerow]['anti'].append(antibio)
            for plusrow in plusdict[genomerow][generow]:
                row += ',' + plusrow

    uniquedict = {}
    for genome in genedict:
        for antibiotic in genedict[genome]['anti']:
            if antibiotic not in uniquedict:
                uniquedict[antibiotic] = []
            for compargenome in genedict:
                if antibiotic not in genedict[compargenome]['anti']:
                    uniquedict[antibiotic].append(genome)
    for antibiotic in uniquedict:
        if len(antibiotic) == 1:
            print antibiotic, uniquedict[antibiotic]
            for genome in uniquedict[antibiotic]:
                genedict[genome]['unique'].append(antibiotic)
    for genomerow in genedict:
        antirow += "%s,\"" % (genomerow)
        for rgene in genedict[genomerow]['gene']:
            antirow += "%s\n" % (rgene)
        antirow = antirow.rstrip()
        antirow += "\",\""
        for ranti in genedict[genomerow]['anti']:
            antirow += "%s\n" % (ranti)
        antirow = antirow.rstrip()
        if len(genedict[genomerow]['unique']) > 0:
            antirow += "\",\""
            for uanti in genedict[genomerow]['unique']:
                antirow += "%s\n" % (uanti)
            antirow = antirow.rstrip()
        antirow += "\"\n"

    with open("%sGeneSeekr_results_%s.csv" % (out, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
        csvfile.write(csvheader)
        csvfile.write(row)
    with open("%sARMI_results_%s.csv" % (out, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
        csvfile.write(antihead)
        csvfile.write(antirow)





parser = argparse.ArgumentParser(description='Antibiotic Resistance Marker Identifier:\n'
                                             'Use to find markers for any bacterial genome')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
parser.add_argument('-i', '--input', required=True, help='Specify input fasta folder')
parser.add_argument('-m', '--marker', required=True, help='Specify antibiotic markers folder')
parser.add_argument('-o', '--output', required=True, help='Specify output folder for csv')
parser.add_argument('-t', '--tab', type=str, required=True, help='tables file location')
args = vars(parser.parse_args())
blaster(args['marker'], args['input'], args['output'])
