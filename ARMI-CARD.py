__author__ = 'mikeknowles'
import json

jsonfile = '/nas/plusdict.json'
json2 = "/nas/Pipeline_development/AntimicrobialResistance/ARMI/tabs/aro3.json"
plusdict = json.load(open(jsonfile))
antidict = json.load(open(json2))

class card():
    '''
    CARD requires a gene # as an input class has three functions:
    card(gene)
    resist will :return and list of antibiotics for a revelant gene
        Includes functionality to trace the depenedencies of a gene complex
        Utilizes recurssion to achieve all possible antibiotics
        .resist(genome)
    anti will :return a list of antibiotics for a given gene
        .anti()
    sens will :return a list of sensitivities for a given gene
        .sens()
    '''
    def __init__(self, index):  # Initialize class and inherit self
        self.index = index  # Defines a gene number that can be passed to functions within the class

    def resist(self, genome=None):  # Begin resist function and import initialized self
        resistlist = []  # Initialize list
        if "resist" in antidict[self.index]:  # If the key "resist" in gene
            if "complex" in antidict[self.index] and genome is not None:
                '''If the key complex in gene defines their are depenedenies'''
                count = 0  # for each resistance set count at zero
                for complex in antidict[self.index]["complex"]:  # Allow for multiple dependencies
                    try:  # some key error exist from card database
                        if plusdict[genome][complex] == ['+']:  # check if dependencies are satisfied
                            count += 1
                    except KeyError:
                        if complex == 3000105:
                            count += 1
                    if len(antidict[self.index]["complex"]) >= count:  # if complexes are satisfied
                        resistlist.extend(antidict[self.index]["resist"])  # extend the list
            else:  # if no complex then just return the list
                resistlist.extend(antidict[self.index]["resist"])
        if "isa" in antidict[self.index]:  # Recursion for parent antibiotic traits
            for depend in antidict[self.index]["isa"]:
                resistlist.extend(card(depend).resist(genome))  # Call self to recurse through the same class
        return resistlist  # return the extended list

    def anti(self):
        if "resist" in antidict[self.index]:
            return antidict[self.index]["resist"]

    def sens(self):
         if "sensitivity" in antidict[self.index]:
            return antidict[self.index]["sensitivity"]

    def function(self):
        if "function" in antidict[self.index]:
            return antidict[self.index]["function"]


class dictbuild():
    '''
    Simple class to build a list or dictionary without repeats
    '''
    def __init__(self, index):
        self.key = index

    def add(self, lst):
        if not self.key:
            self.key = lst
        else:
            for drug in lst:
                if drug not in self.key:
                    self.key.append(drug)
        # self.key = sorted(self.key, key=lambda s: s.lower())
        return self.key



outputdict = {}
for genome in sorted(plusdict):  # iterate through plus dict
    if "/sig/EcoliAll/OLC714.fas" not in outputdict or "/sig/EcoliAll/OLC714" not in genome:
        outputdict[genome] = {"resist": [], "sensitivity": [], "genes": []}
        resistance = outputdict[genome]["resist"]
        for gene in plusdict[genome]:
            # refgene = plusdict[genome][gene]
            if plusdict[genome][gene]:
                resist = card(gene).resist(genome)  # check resistances
                sens = card(gene).sens()  # check sensitivities
                if resist is not None:
                   outputdict[genome]["resist"] = dictbuild(outputdict[genome]["resist"]).add(resist)
                   outputdict[genome]["genes"].append(gene)
                   # print json.dumps(outputdict[genome], sort_keys=True, indent=4, separators=(',', ': '))
                if sens is not None:
                    outputdict[genome]["sensitivity"] = dictbuild(outputdict[genome]["sensitivity"]).add(sens)
json.dump(outputdict, open('/nas/knowlesm/test2.json', 'w'), sort_keys=True, indent=4, separators=(',', ': '))
antilist = []


for gene in antidict:  # build hearder list
    resistances = card(gene).anti()
    if resistances is not None:
        for resist in resistances:
            if resist not in antilist:
                antilist.append(resist)
antihead = "Genome"
drugcounter = {}
antilist = sorted(antilist, key=lambda s: s.lower())  # sort header case insensitive
for anti in antilist:
    antihead += ",\"%s\"" % anti
    drugcounter[anti] = 0

olc = 0
antistr = ""

''' Build csv string '''
for genome in sorted(outputdict):
    genomename = "\n\"%s\"" % genome.split("/")[-1].split(".")[0]
    antistr += genomename
    genomecount = 0
    for drug in antilist:
        if drug in outputdict[genome]["resist"]:
            antistr += ",+"
            drugcounter[drug] += 1
            genomecount += 1
        else:
            antistr += ",-"
    antistr += ",%i" % genomecount
antihead += "\nCount"
for drug in antilist:
    antihead += ",%i" % drugcounter[drug]
antihead += antistr



open('/nas/knowlesm/test2.csv', 'w').write(antihead)
