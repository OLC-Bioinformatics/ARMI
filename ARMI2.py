from collections import defaultdict
import os, time, argparse, json, GeneSeekrv2


def defdict(file):
    ndict = defaultdict(list)
    with open(args['tab']+file) as tab:
        for line in tab:
            tline = line.strip().split(None, 1)
            ndict[tline[0]].append(tline[1])
    return ndict

def blaster(path, targets, out, threshold):
    jsonfile = '%splusdict.json' % targets
    print jsonfile
    if os.path.isfile(jsonfile):
        plusdict = json.load(open(jsonfile))

    else:
        plusdict = GeneSeekrv2.blaster(path, targets, out, 'ARMI2')
        json.dump(plusdict, open(jsonfile, 'w'), sort_keys=True, indent=4, separators=(',', ': '))
    count = 80
    require = defdict("require.tab")
    antidict = defdict("typeResis.tab")
    genedict = defaultdict(dict)
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
                            if plusdict[genomerow][propercasegene] != 0:
                                tempcount += 1
            elif plusdict[genomerow][generow] != 0:
                # print plusdict[genomerow][generow]
                tempcount = 1
            else:
                tempcount = 0
                # print tempcount, len(require[generow.lower()]), require[generow.lower()]
            if tempcount == len(require[generow.lower()]) | 1 and plusdict[genomerow][generow] != 0:
                genedict[genomerow]['gene'].append(generow)
                # print genomerow, generow, plusdict[genomerow][generow]
                for antibio in antidict[generow.lower()]:
                    if antibio not in genedict[genomerow]['anti']:
                        genedict[genomerow]['anti'].append(antibio)
            # for plusrow in plusdict[genomerow][generow]:
            row += ',' + str(plusdict[genomerow][generow])

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
    json.dump(uniquedict, open('%s/../unique.json' % path, 'w'), sort_keys=True, indent=4, separators=(',', ': '))
    json.dump(antidict, open('%s/../anti.json' % path, 'w'), sort_keys=True, indent=4, separators=(',', ': '))
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