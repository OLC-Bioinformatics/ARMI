__author__ = 'mike'

import json
import itertools
from multiprocessing import Process, Manager, Pool

test = json.load(open("/Users/mike/Google Drive/CFIA/ARMI/test2.json"))

amrlist = []
combination = []

with open("/Users/mike/Google Drive/CFIA/ARMI/test2.csv") as csv:
    antilist = csv.readline().replace('\"','').rstrip().split(",")[2:]
    countlst = csv.readline().replace('\"','').rstrip().split(",")[2:]
    for i in range(len(countlst)):
        if int(countlst[i]) <= 215:
            amrlist.append(antilist[i])
    print amrlist
    print countlst
    print len(amrlist), len(antilist)
    combination = list(itertools.combinations(antilist,2))
    print len(combination)

'''
Use combinations in targeted way to retrieve more unique results


'''
# def chkcock(cockdict, somelist):
genome = "/sig/EcoliAll/OLC714.fas"
print test[genome]["resist"]
    # for genome in test:
print genome
for cocktail in combination:
    print
    if {cocktail}.issuperset(test[genome]["resist"]):
        if cocktail not in cockdict:
            cockdict[cocktail] = []
        cockdict[cocktail].append(genome)

# if __name__ == '__main__':
#     manager = Manager()
#     d = manager.dict()
#     l = manager.list(range(10))
#
#     p = Process(target=chkcock, args=(d, l))
#     p.start()
#     p.join()
#     print l
#     print d
