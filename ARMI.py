__author__ = 'mikeknowles'
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from subprocess import call
import os
import glob
import time
import sys
import threading
import Queue
import shlex
import argparse

queue = Queue.Queue()

class makeblastdb(threading.Thread):
    """Threaded makeblastdb"""
    def __init__(self, queue):
      threading.Thread.__init__(self)
      self.queue = queue

    def run(self):
        while True:
            #grabs fastapath from queue
            fastapath = self.queue.get()
            #makeblastdb if not exist
            if not os.path.isfile(fastapath + ".nhr"):
                call(shlex.split("makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, fastapath)))
            else:
                pass
            # signals to queue job is done
            self.queue.task_done()

def makedbthreads(path):
    ''' Setup and create threads for class'''
    fastas = glob.glob(path + "*.fasta")
    for i in range(len(fastas)):
        sys.stdout.write('.')
        t = makeblastdb(queue)
        t.setDaemon(True)
        t.start()
    for fasta in fastas:
        queue.put(fasta)
    #wait on the queue until everything has been processed
    queue.join()

def blaster(path, out):
    print "[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S"))
    makedbthreads(path)
    print "[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S"))


parser = argparse.ArgumentParser(description='Antibiotic Resistance Marker Identifier:\n'
                                             'Use to find markers for any bacterial genome')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
parser.add_argument('-i', '--input', required=True, help='Specify input fasta folder')
parser.add_argument('-m', '--marker', required=True, help='Specify antibiotic markers folder')
parser.add_argument('-o', '--output', required=True, help='Specify output folder for csv')
args = vars(parser.parse_args())
blaster(args['input'], args['output'])
