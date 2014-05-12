__author__ = 'mike knowles'
#!/usr/bin/python
from Bio import Entrez, SeqIO
import os
import time


def FileOpen(FileLocation):
    print "File location is: %s"  % (FileLocation)
    tab = open(os.path.expanduser(FileLocation), 'r')
    for line in tab:
        cell = line.split("\t")
        current = time.strftime("%H:%M:%S")
        # cell[0]  is the accession number of a bacterial genome in NCBI.
        # cell[6]  is the query sequence length.
        # cell[10] is the start coordinate of query sequence in the HSP.
        # cell[11] is the end coordinate of query sequence in the HSP.
        print "[%s] For the query genome: %s with the length %s has the coordinates: %s-%s"\
              % (current, cell[0], cell[6], cell[10], cell[11])
        time.sleep(5)
        GenomeLookup(cell[0], cell[10], cell[11])
    tab.close()
    return

def GenomeLookup(ac,start, stop):
    Entrez.email = "email@example.com"
    search = Entrez.esearch(db="nuccore",
                            term=str(ac),
                            )
    result = Entrez.read(search)
    handle = Entrez.efetch(db="nuccore",
                           id=result["IdList"][0],
                           rettype="fasta",
                           strand=1,
                           seq_start=start,
                           seq_stop=stop)
    record = SeqIO.read(handle, "fasta")
    handle.close()
    current = time.strftime("%H:%M:%S")
    print "[%s] id is %s" % (current, result["IdList"][0])
    print "[%s] %s" % (current, record.seq)
filehandle = raw_input("Enter the table file: ") + "/genomeblast.tab"
FileOpen(filehandle)
