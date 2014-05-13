__author__ = 'mike knowles'
from Bio import Entrez, SeqIO
import os
import time
import re


def FileOpen(FileLocation):
    current = time.strftime("%H:%M:%S")
    prodict = {}
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
            # time.sleep(5)
            # GenomeLookup(aclength.group(1), aclength.group(2), aclength.group(3), progene.group(2))
            if aclength.group(1) not in prodict:
                prodict[aclength.group(1)] = {progene.group(1): {'start': aclength.group(2), 'stop': aclength.group(3)}}
            else:
                prodict[aclength.group(1)][progene.group(1)] = {'start': aclength.group(2), 'stop': aclength.group(3)}
                print prodict[aclength.group(1)][progene.group(1)].items()
            # print "[%s] For the protein %s & resistance gene %s belonging to query genome: " \
            #       "%s with the length %s has the coordinates: %s-%s" \
            #       % (current, progene.group(1), progene.group(2), aclength.group(1),
            #      int(aclength.group(3))-int(aclength.group(2)), aclength.group(2), aclength.group(3))
        elif progene:
            # print "[%s] No genome for %s protein & %s resistance gene" % (current, progene.group(1), progene.group(2))
            continue
        else:
            print "[%s] Failure" % (current)
        if len(prodict[aclength.group(1)].keys()) < 1:
            print prodict[aclength.group(1)].keys()
    tab.close()
    # dictprinter(prodict)  # Test dictionary
    print (len(prodict.items()))
    return

def dictprinter(obj):
    for x in obj:
        print ("\n" + x)
        for y in obj[x]:
            print (y)
            for z in obj[x][y]:
                print (z + ":" + obj[x][y][z])

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