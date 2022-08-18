"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file = open(filename,"r")
    text = file.read()
    string = text.replace("\n","")    
    return string

'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    dna = dna.replace("T","U")
    dna = dna[startIndex:]
    codon =  ""
    ends = ["UAA", "UAG", "UGA"]
    codons = []
    for i in dna:
        codon += i
        if len(codon)==3:
            codons.append(codon)
            if codon in ends:
                break
            codon = ""
    return codons
    


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f = open(filename)
    data = json.load(f)  
    codon = {}
    for i,j in data.items():        
        for k in j:
            if "T" in k: 
                k = k.replace("T","U")
            codon[k]=i        
    return codon


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD): 
    protein =[]    
    for i in codons:                   
        protein.append(codonD[i]) 
        if protein[0]=="Met":
            protein[0] = "Start" 
    return protein


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna =readFile(dnaFilename)
    codon = makeCodonDictionary(codonFilename)    
    synthesize = []
    index = 0 
    un_used_bases = 0
    while index<len(dna):
        codon1 = dna[index:index+3]
        if codon1=="ATG":
            rna = dnaToRna(dna,index)
            protein = generateProtein(rna,codon)
            synthesize.append(protein)
            index += 3*len(protein)
        else:
            index += 1
            un_used_bases += 1  
    print("total number of bases",len(dna),"unused-base count",un_used_bases,"total number of proteins synthesized",len(synthesize),"\n")                 
    return synthesize

def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    common = [i for i in proteinList1 if i in proteinList2]
    
    return common


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    
    return sum(proteinList,[])


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    amino = {i:aaList.count(i) for i in aaList }
    return amino


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):    
    combine1 = combineProteins(proteinList1)
    combine2 =  combineProteins(proteinList2)
    amino1 = aminoAcidDictionary(combine1)
    amino2 = aminoAcidDictionary(combine2)    
    unique = list(set(combine1+combine2)) 
    freq = []   
    for k in unique:
        f1=0
        f2=0        
        if k!="Start" and k!="Stop":
            if k in amino1:
                f1 = amino1[k]/len(combine1)
            else:
                f1 = 0
            if k in amino2:
                f2 = amino2[k]/len(combine2)
            else:
                f2 = 0
                
            d = abs(f1-f2)       
            if d>cutoff:
                freq.append([k,round(f1,4),round(f2,4)])               
    return freq


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    list1=[]
    for i in commonalities:
        i.remove("Start")
        i.remove("Stop")
        if len(i)>1:
            new='-'.join(i)
            list1.append([new])
        else:
            if i not in list1:
                list1.append(i)
    new=sorted(list1)
    print("The following proteins occurred in both DNA Sequences:","\n")
    for i in new:
        for j in i:
            print(j)
    print("\n")
    print("The following amino acids occurred at very different rates in the two DNA sequences:","\n")
    for i in differences:
        print(i[0]+":"+str(round(i[1]*100,2))+"%"+ " in seq1 ,"+str(round(i[2]*100,2))+"%"+"in seq2")
    print("\n")
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    combine1 = combineProteins(proteinList1)
    combine2 = combineProteins(proteinList2)
    amino = sorted(set(combine1+combine2))    
    return amino



'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    amino = combineProteins(proteinList)
    amino_d = aminoAcidDictionary(amino)
    freq = [amino_d[i]/len(amino) if i in amino else 0  for i in labels ]
    
        
    return freq


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()

    ## Uncomment these for Week 2 ##
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    

    ## Uncomment these for Week 3 ##
   
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
