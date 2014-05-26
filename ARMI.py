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
import subprocess, os, glob, time, sys, shlex, argparse, re

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
    for i in range(len(fastas) * len(genomes)):
        threads = Thread(target=runblast, args=(blastqueue,))
        threads.setDaemon(True)
        threads.start()
    blastpath = []
    for genome in genomes:
        for fasta in fastas:
            gene = re.search('\/(\w+)\.fasta', fasta)
            path = re.search('(.+)\/(.+?)\.fasta', genome)
            out = "%s/../tmp/%s.%s.xml" % (path.group(1), path.group(2), gene.group(1))
            blastpath.append(out)
            if not os.path.isfile(out):
                blastqueue.put((genome, fasta, out))
                dotter()
    #wait on the queue until everything has been processed
    blastqueue.join()
    return blastpath


def blastparser(parsequeue):
    while True:
        hsp = parsequeue.get()
        if hsp.identities == hsp.align_length:
            plus = '+'
        elif hsp.align_length > hsp.identities >= (hsp.align_length * 0.9):
            print "%s %s" % (hsp.align_length, hsp.identities)
            perid = int((float(hsp.identities) / hsp.align_length) * 100)
            print(perid)
            plus = '(%s%%)' % (perid)
        else:
            plus = '-'
        parsequeue.task_done()
        parsequeue.put(plus)
        return plus


def parsethreader(xml):
    dotter()
    numhsp = 0
    with open(xml) as file:
        numhsp = sum(line.count('<Hsp>') for line in file)
    if numhsp != 0:
        handle = open(xml)
        records = NCBIXML.parse(handle)
        for i in range(numhsp):
            threads = Thread(target=blastparser, args=(parsequeue,))
            threads.setDaemon(True)
            threads.start()
        for record in records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    parsequeue.put(hsp)
                    plus = parsequeue.get()
                    print plus


def blaster(path, targets, out):
    fastas = glob.glob(path + "*.fasta")
    genomes = glob.glob(targets + "*.fasta")
    sys.stdout.write("[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S")))
    makedbthreads(fastas)
    print "\n[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S"))
    print "[%s] Now performing BLAST database searches" % (time.strftime("%H:%M:%S"))
    blastpath = blastnthreads(fastas, genomes)
    for xml in blastpath:
     parsethreader(xml)
    parsequeue.join()



parser = argparse.ArgumentParser(description='Antibiotic Resistance Marker Identifier:\n'
                                             'Use to find markers for any bacterial genome')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
parser.add_argument('-i', '--input', required=True, help='Specify input fasta folder')
parser.add_argument('-m', '--marker', required=True, help='Specify antibiotic markers folder')
parser.add_argument('-o', '--output', required=True, help='Specify output folder for csv')
args = vars(parser.parse_args())
blaster(args['marker'], args['input'], args['output'])
