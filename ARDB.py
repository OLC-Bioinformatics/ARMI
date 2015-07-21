__author__ = 'mike'
import json


antidict = json.load(open("/Users/mike/Google Drive/CFIA/tabs/aro3.json"))

# antidict = {}

def addgene(gene, dict):
    dict[gene] = {'resist': [], 'name': gene, 'function': '', 'complex': [], 'db': "ardb"}


for line in open("/Users/mike/Google Drive/CFIA/tabs/typeResis.tab"):
    linelist = line.split("\t")
    if linelist[0] not in antidict:
        addgene(linelist[0], antidict)

    antidict[linelist[0]]['resist'].append(linelist[1].rstrip())


for line in open("/Users/mike/Google Drive/CFIA/tabs/require.tab"):
    linelist = line.split("\t")
    if linelist[0] not in antidict:
        antidict[linelist[0]] = {'resist': [], 'name': linelist[0], 'function': [], 'complex': []}
    antidict[linelist[0]]['complex'].append(linelist[1].rstrip())

for line in open("/Users/mike/Google Drive/CFIA/tabs/classinfo.tab"):
    linelist = line.split("\t")
    if linelist[0] not in antidict:
        antidict[linelist[0]] = {'resist': [], 'name': linelist[0], 'function': [], 'complex': []}
    antidict[linelist[0]]['function'] = linelist[1].rstrip()

json.dump(antidict, open("aro.json", "w"), sort_keys=True, indent=4, separators=(',', ': '))
