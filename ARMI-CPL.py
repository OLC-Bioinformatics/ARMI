__author__ = 'mikeknowles'
'''Iterate through dictionary to find unique carbon and antibiotics'''
from collections import defaultdict
import json, textwrap
aro = ""
cardict = {}
for line in open("/nas/Pipeline_development/AntimicrobialResistance/ARMI/tabs/aro.obo"):
    # print line, "[Term]" in line
    if "id: ARO" == line[:7]:
        aro = line[8:].rstrip()
        cardict[aro] = defaultdict(list)
        # cardict[aro] = {"name": {}, "resist": [], "complex": [], "function": {}, "sensitivity": []}
        # cardict[aro]["name"] = {}
    if "name: " == line[:6]:
        cardict[aro]["name"] = line[6:].rstrip()
    if "is_a: ARO:" in line:
        cardict[aro]["isa"].append(line[10:18].rstrip())
        cardict[aro]["function"] = line[20:].rstrip()
    if "relationship: confers_resistance_to_drug" in line:
        cardict[aro]["resist"].append(line[55:].rstrip())
    elif "relationship: confers_resistance_to " in line:
        cardict[aro]["resist"].append(line[50:].rstrip())
    if "relationship: part_of" in line:
        cardict[aro]["complex"].append(line[26:33].rstrip())
    if "relationship: targeted_by_drug" in line:
        cardict[aro]["sensitivity"].append(line[45:].rstrip())
# deletelist = []
# for record in cardict:
#     if "resist" not in cardict[record]:
#         deletelist.append(record)
# for card in deletelist:
#     del cardict[card]
json.dump(cardict,open("/nas/Pipeline_development/AntimicrobialResistance/ARMI/tabs/aro3.json",'w')
          , indent=4, separators=(',', ': '))



# genelist = {}
# for line in open("/nas/Pipeline_development/AntimicrobialResistance/ARMI/tabs/AROtags.txt"):
#     linelist = line.rstrip().split("\t")
#     splitgene = linelist[0].split(".")
#     gene = ".".join(splitgene[:-1])
#     # print gene
#     genelist[gene] = linelist[2][4:]
# print genelist
# handle = open("/dev/null", 'w')
# for line in open("/nas/Pipeline_development/AntimicrobialResistance/ARMI/tabs/AR-genes.fa"):
#     if ">" == line[0]:
#         try:
#             aro = genelist[line.split(" ")[0][1:]]
#         except KeyError:
#             # print line
#             aro = line.split("ARO:")[-1][:7]
#             if aro == "1000001":
#                 aro = line.split("ARO:")[1][:7]
#             print aro
#         handle.close()
#         handle = open("/nas/Pipeline_development/"
#                       "AntimicrobialResistance/ARMI/card/%s.fasta" % aro, 'a')
#         handle.write(line)
#     else:
#         handle.write("\n".join(textwrap.wrap(line, width=70)) + "\n")