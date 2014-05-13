__author__ = 'mike knowles'
from Bio import Entrez, SeqIO
import os
import time
import re


def FileOpen(FileLocation):
    current = time.strftime("%H:%M:%S")
    print "[%s] File location is: %s" % (current, FileLocation)
    tab = open(os.path.expanduser(FileLocation), 'r')
    for line in tab:
        progene = re.search('^(\S+)\t(\S+)', line)
        aclength = re.search('\tC?_?(\S+?)\|<?(\d+)\|>?(\d+)+', line)
        current = time.strftime("%H:%M:%S")
        # progene.group(1) is the accession number of a bacterial protein in NCBI.
        # progene.group(2) is the name of the resistance gene
        # aclength.group(1)  is the accession number of a bacterial genome in NCBI.
        # Calculate difference for the query sequence length.
        # aclength.group(2) is the start coordinate of query sequence in the HSP.
        # aclength.group(3) is the end coordinate of query sequence in the HSP.
        if aclength:
            print "[%s] For the protein %s & resistance gene %s belonging to query genome: " \
                  "%s with the length %s has the coordinates: %s-%s" \
                  % (current, progene.group(1), progene.group(2), aclength.group(1),
                     int(aclength.group(3))-int(aclength.group(2)), aclength.group(2), aclength.group(3))
            # time.sleep(5)
            # GenomeLookup(aclength.group(1), aclength.group(2), aclength.group(3), progene.group(2))
        elif progene:
            print "[%s] No genome for %s protein & %s resistance gene" % (current, progene.group(1), progene.group(2))
            prodict = {progene.group(1): {aclength.group(1): {'start': aclength.group(2), 'stop': aclength.group(3)}}}
        else:
            print "[%s] Failure" % (current)
    tab.close()
    dictprinter(prodict)  # Test dictionary
    return

def dictprinter(obj):
    for key, obj in obj.items():
    print(key)
    for attribute, value in obj.items():
        print('{} : {}'.format(attribute, value))

def GenomeLookup(ac,start, stop, gene):
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
filehandle = "/nas/Pipeline_development/AntimicrobialResistance/ardbAnno1.0/tabs/abrg.tab" # raw_input("Enter the table files folder: ") + "abrg.tab"
FileOpen(filehandle)