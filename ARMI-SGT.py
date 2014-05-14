__author__ = 'mike knowles'
from Bio import Entrez, SeqIO
import os
import time
import re
import sys
import argparse

def FileOpen(FileLocation):
    current = time.strftime("%H:%M:%S")
    prodict = {}
    print "[%s] File location is: %s" % (current, FileLocation)
    tab = open(os.path.expanduser(FileLocation), 'r')
    for line in tab:
        progene = re.search('^(\S+)\t(\S+)', line)
        aclength = re.search('\tC?_?(\S+?)\|<?(\d+)\|>?(\d+)$', line)
        current = time.strftime("%H:%M:%S")
        # progene.group(1) is the accession number of a bacterial protein in NCBI.
        # progene.group(2) is the name of the resistance gene
        # aclength.group(1)  is the accession number of a bacterial genome in NCBI.
        # Calculate difference for the query sequence length.
        # aclength.group(2) is the start coordinate of query sequence in the HSP.
        # aclength.group(3) is the end coordinate of query sequence in the HSP.

        if aclength:
            try:
                if aclength.group(2) in prodict[aclength.group(1)]:
                    prodict[aclength.group(1)][aclength.group(2)][aclength.group(3)] = {progene.group(2)}
                else:
                    prodict[aclength.group(1)][aclength.group(2)] = {aclength.group(3): {progene.group(2)}}
            except:
                prodict[aclength.group(1)] = {aclength.group(2): {aclength.group(3): {progene.group(2)}}}
                # print prodict[aclength.group(1)][progene.group(1)].items()
            # print "[%s] For the protein %s & resistance gene %s belonging to query genome: " \
            #       "%s with the length %s has the coordinates: %s-%s" \
            #       % (current, progene.group(1), progene.group(2), aclength.group(1),
            #      int(aclength.group(3))-int(aclength.group(2)), aclength.group(2), aclength.group(3))
        else:
            pass
            # print "[%s] Failure %s" % (current, progene.group(1))
    tab.close()
    # dictprinter(prodict)  # Test dictionary
    decider(prodict)
    print (len(prodict.items()))
    return

def dictprinter(obj):
    for x in obj:
        print ("\n" + x)
        for y in obj[x]:
            print (y)
            for z in obj[x][y]:
                print (z + ":" + obj[x][y][z])

def decider(prodict):
    '''
    Decide if full genome or genome fragment should be downloaded with GenomeLookup
    try:
        abrgfasta = sys.argv[1]
    except:
        abrgfasta = '''"/nas/Pipeline_development/AntimicrobialResistance/ARMI/db/abrg.fa" ''' raw_input("Enter the save location:")
    '''
    abrg = open(args['output'], 'w')
    for genomeac in prodict:
        current = time.strftime("%H:%M:%S")
        try:
            if counter > 0:
                print "[%s] Just finished writing %s sequences" %(current, counter)
        except UnboundLocalError:
            pass
        if len(prodict[genomeac].items()) == 1:
            single = True
            # print len(prodict[genomeac].items())
        else:
            # print (prodict[genomeac].items())
            counter = 0
            single = False
            record = GenomeLookup(genomeac)
        for begin in prodict[genomeac]:
            for end in prodict[genomeac][begin]:
                for gene in prodict[genomeac][begin][end]:
                    try:
                        if single:
                            counter = 0
                            record = GenomeLookup(genomeac, begin, end)
                            # print ">%s (%s)\n%s" % (record.description, gene, record.seq)
                            abrg.write(">%s (%s)\n%s\n" % (record.description, gene, record.seq))
                        else:
                            abrg.write(">%s (%s) %s-%s\n%s\n" % (record.description, gene, begin, end, record.seq[int(begin):int(end)]))
                            counter += 1
                    except AttributeError:
                        print "[%s] Unable to print %s to file" % (current, genomeac)
    abrg.close()




def GenomeLookup(ac, start=False, stop=False):
    time.sleep(1)
    Entrez.email = args['email']
    current = time.strftime("%H:%M:%S")
    search = Entrez.esearch(db="nuccore",
                            term=str(ac),
                            )
    result = Entrez.read(search)
    try:
        if start == False:
            handle = Entrez.efetch(db="nuccore",
                           id=result["IdList"][0],
                           rettype="fasta")
            record = SeqIO.read(handle, "fasta")
            handle.close()
        else:
            handle = Entrez.efetch(db="nuccore",
                                   id=result["IdList"][0],
                                   rettype="fasta",
                                   strand=1,
                                   seq_start=start,
                                   seq_stop=stop)
            record = SeqIO.read(handle, "fasta")
            handle.close()
        print "[%s] ID is %s for accession %s" % (current, result["IdList"][0], ac)
        return record
    except:
        print "[%s] Unable to locate genome %s, Unexpected error: %s" % (current, ac, sys.exc_info()[0])

"Parser for argurments"
parser = argparse.ArgumentParser(description='Convert table file to fasta genomes given coordinates')
parser.add_argument('--version', action='version', version='%(prog)s v0.3')
parser.add_argument('-t', '--tab', type=str, required=True, help='arbg.tab file location')
parser.add_argument('-e', '--email', required=True, help='A valid email address is required')
parser.add_argument('-o', '--output', required=True, help='Specify output file')
args = vars(parser.parse_args())
FileOpen(args['tab'])