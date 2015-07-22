__author__ = 'mikeknowles'
from glob import glob
from re import match
from ARMICARD import decipher
import os, time, argparse, json, GeneSeekrv2




def blaster(path, targets, out, threshold, db):
    if db == "both":
        db = ['ardb', 'card']
    else:
        db = [db]
    jsonfile = '%splusdict.json' % targets
    print jsonfile
    if os.path.isfile(jsonfile):
        plusdict = json.load(open(jsonfile))

    else:
        markers = glob(path + "/*")
        for marker in markers:
            cardcheck = match("^\d{7}$", marker)
            if db == 'ardb' and cardcheck is not None:
                markers.remove(marker)
            elif db == 'card' and cardcheck is None:
                markers.remove(marker)

        plusdict = GeneSeekrv2.blaster(markers, targets, out, 'ARMI2')
        json.dump(plusdict, open(jsonfile, 'w'), sort_keys=True, indent=4, separators=(',', ': '))
    antidict = json.load("aro.json")
    decipher(plusdict, antidict)



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
parser.add_argument('-d', '--db', default='both', help='Specify antibiotic markers database')
args = vars(parser.parse_args())
blaster(args['marker'], args['input'], args['output'], args['cutoff'], args['db'])