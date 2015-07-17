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
import subprocess, os, glob, time, sys, shlex, argparse, re, threading, json, mmap

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
            # print "makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, fastapath)
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
    '''Setup and create  threads for blastn and xml path'''
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
            blastpath[out] = {path.group(2): (gene.group(1),)}
            if not os.path.isfile(out):
                blastqueue.put((genome, fasta, out))
                dotter()
        #wait on the queue until everything has been processed
        blastqueue.join()
    return blastpath

plusdict = {}


def blastparser(alignment, hsp, genomes):
    global plusdict
    if hsp.identities == alignment.length:
        plus = '+'
    elif alignment.length > hsp.identities >= (alignment.length * 0.7):
        perid = int((float(hsp.identities) / alignment.length) * 100)
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
        handle = open(xml, 'r')
        mm = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
        handle.close()
        records = NCBIXML.parse(mm)
        for record in records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    blastparser(alignment, hsp, genomes)
        mm.close()
    else:
        for genome in genomes:
            for gene in genomes[genome]:
                try:
                    if gene not in plusdict[genome]:
                        plusdict[genome][gene] = {'-'}
                except:
                    plusdict[genome] = {gene: {'-'}}


def parsethreaderr(xml, genomes):
    global plusdict
    dotter()
    numhsp = 0
    with open(xml) as file:
        numhsp = sum(line.count('<Hsp>') for line in file)
    if numhsp != 0:
        handle = open(xml, 'r')
        records = NCBIXML.parse(handle)
        for record in records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    blastparser(alignment, hsp, genomes)
        handle.close()
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

def blaster(path, targets, out, threshold):
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
    if os.path.isfile('%sblastxmldict.json' % targets):
        print "[%s] Loading BLAST data from file" % (time.strftime("%H:%M:%S"))
        blastpath = json.load(open('%sblastxmldict.json' % targets))
    else:
        print "[%s] Now performing BLAST database searches" % (time.strftime("%H:%M:%S"))
        # make blastn threads and retrieve xml file locations
        blastpath = blastnthreads(fastas, genomes)
        json.dump(blastpath, open('%sblastxmldict.json' % targets, 'w'), sort_keys=True, indent=4, separators=(',', ': '))
    count = 80
    require = defdict("require.tab")
    antidict = defdict("typeResis.tab")
    genedict = defaultdict(dict)
    start = time.clock()
    for xml in blastpath:
        parsethreaderr(xml, blastpath[xml])
    io = time.clock() - start
    start = time.clock()
    for xml in blastpath:
        parsethreader(xml, blastpath[xml])
    mmaped = time.clock() - start
    print "\n[%s] Runtime for mmap method: %ss I/O method: %ss" % (time.strftime("%H:%M:%S"), mmaped, io)
    csvheader = 'Strain'
    row = ""
    antirow = ""
    uniquerow = ""
    rowcount = 0
    antihead = "Strain,Gene,\"Resistance Antibiotic\",\"Unique Resistance\"\n"
    uniquehead = "Strain"
    for genomerow in plusdict:
        genedict[genomerow] = {'anti': [], 'gene': [], 'unique': []}
        row += "\n" + genomerow
        rowcount += 1
        for generow in sorted(plusdict[genomerow]):
            if rowcount <= 1:
                csvheader += ',' + generow
            tempcount = 1

            if generow.lower() in require:
                for othergene in require[generow.lower()]:
                    for propercasegene in plusdict[genomerow]:
                        if othergene == propercasegene.lower():
                            # print propercasegene, othergene, require[generow.lower()], tempcount, len(require[generow.lower()]), plusdict[genomerow][generow]
                            if plusdict[genomerow][propercasegene] != {'-'}:
                                tempcount += 1
            elif plusdict[genomerow][generow] != {'-'}:
                # print plusdict[genomerow][generow]
                tempcount = 1
            else:
                tempcount = 0
                # print tempcount, len(require[generow.lower()]), require[generow.lower()]
            if tempcount == len(require[generow.lower()]) | 1 and plusdict[genomerow][generow] != {'-'}:
                genedict[genomerow]['gene'].append(generow)
                # print genomerow, generow, plusdict[genomerow][generow]
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
                if antibiotic not in genedict[compargenome]['anti'] and genome not in uniquedict[antibiotic]:
                    uniquedict[antibiotic].append(genome)
    for antibiotic in sorted(uniquedict):
        uniquehead += "," + antibiotic
        # print len(uniquedict[antibiotic]), uniquedict[antibiotic]
        if len(uniquedict[antibiotic]) <= threshold:
            for genome in uniquedict[antibiotic]:
                if len(uniquedict[antibiotic]) != 1:
                    antithres = "%s (%s)" % (antibiotic, len(uniquedict[antibiotic]))
                else:
                    antithres = antibiotic
                genedict[genome]['unique'].append(antithres)
    uniquehead += "\n"
    for genomerow in genedict:
        antirow += "%s,\"" % (genomerow)
        uniquerow += "%s" % (genomerow)
        for rgene in genedict[genomerow]['gene']:
            antirow += "%s\n" % (rgene)
        antirow = antirow.rstrip()
        antirow += "\",\""
        for ranti in genedict[genomerow]['anti']:
            antirow += "%s\n" % (ranti)
        for antibiotic in sorted(uniquedict):
            if antibiotic in genedict[genomerow]['anti']:
                uniquerow += ",+"
            else:
                uniquerow += ",-"
        antirow = antirow.rstrip()
        uniquerow += "\n"
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
    with open("%sARMI_unique_%s.csv" % (out, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
        csvfile.write(uniquehead)
        csvfile.write(uniquerow)




parser = argparse.ArgumentParser(description='Antibiotic Resistance Marker Identifier:\n'
                                             'Use to find markers for any bacterial genome')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
parser.add_argument('-i', '--input', required=True, help='Specify input fasta folder')
parser.add_argument('-m', '--marker', required=True, help='Specify antibiotic markers folder')
parser.add_argument('-o', '--output', required=True, help='Specify output folder for csv')
parser.add_argument('-t', '--tab', type=str, required=True, help='tables file location')
parser.add_argument('-c', '--cutoff', type=int, default=1, help='Threshold for maximum unique bacteria'
                                                                ' for a single antibiotic')
args = vars(parser.parse_args())
blaster(args['marker'], args['input'], args['output'], args['cutoff'])
