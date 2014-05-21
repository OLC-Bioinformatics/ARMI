__author__ = 'mikeknowles'
""" Includes threading found in examples:
http://www.troyfawkes.com/learn-python-multithreading-queues-basics/
http://www.ibm.com/developerworks/aix/library/au-threadingpython/
https://docs.python.org/2/library/threading.html
"""
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from threading import Thread
from Queue import Queue
import subprocess, os, glob, time, sys, shlex, argparse, re


def makeblastdb(queue):
    while True:
        #grabs fastapath from queue
        fastapath = queue.get()
        #makeblastdb if not exist
        nhr = "%s.nhr" % (fastapath)
        if not os.path.isfile(str(nhr)):
            print "makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, fastapath)
            subprocess.Popen(shlex.split("makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, fastapath)),
                             stdout=open(os.devnull, 'wb'))
            sys.stdout.write('.')
        else:
            pass
        # signals to queue job is done
        queue.task_done()

queue = Queue()
queue2 = Queue()

def makedbthreads(fastas):
    ''' Setup and create threads for class'''
    for i in range(len(fastas)):
        threads = Thread(target=makeblastdb, args=(queue,))
        threads.setDaemon(True)
        threads.start()
    for fasta in fastas:
        queue.put(fasta)
    #wait on the queue until everything has been processed
    queue.join()

def runblast(queue2):
    while True:
        genomefasta = queue2.get()
        for key in genomefasta:
            for db in key:
                print db
                gene = re.search('\/(\w+)\.fasta', db)
                blastn = NcbiblastnCommandline(query=key, db=db, evalue=1e-40, out=gene.group(1))
                print blastn
                subprocess.Popen(shlex.split(blastn))
        queue2.task_done()

def blastnthreads(fastas, genomes):
    for i in range(len(fastas) * len(genomes)):
        threads = Thread(target=runblast, args=(queue2,))
        threads.setDaemon(True)
        threads.start()
    genomefasta = {}
    for genome in genomes:
        genomefasta.setdefault(genome, []).append(fastas)
        queue2.put(genomefasta)
    #wait on the queue until everything has been processed
    queue2.join()

def blaster(path, targets, out):
    fastas = glob.glob(path + "*.fasta")
    genomes = glob.glob(targets + "*.fasta")
    sys.stdout.write("[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S")))
    makedbthreads(fastas)
    print "\n[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S"))
    print "[%s] Now performing BLAST database searches" % (time.strftime("%H:%M:%S"))
    blastnthreads(fastas, genomes)




parser = argparse.ArgumentParser(description='Antibiotic Resistance Marker Identifier:\n'
                                             'Use to find markers for any bacterial genome')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
parser.add_argument('-i', '--input', required=True, help='Specify input fasta folder')
parser.add_argument('-m', '--marker', required=True, help='Specify antibiotic markers folder')
parser.add_argument('-o', '--output', required=True, help='Specify output folder for csv')
args = vars(parser.parse_args())
blaster(args['marker'], args['input'], args['output'])
