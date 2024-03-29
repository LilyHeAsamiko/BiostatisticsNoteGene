import pandas as pd
import numpy as np

def RF(File):
    '''Read a txt file and output its content as a string.'''
    f =open(File, "r")
    s= f.read()
    s=s.replace("\n","")
    s=s.replace("\r","")
    return s

def RFAS(File):
    '''Read a fasta file and output its content as three strings:seq,startIndex,stopIndex.'''
    #with open(File, "r")as f:s= f.read()
    f =open(File, "r")
    s= f.read()
    s=s.replace("\n","").replace("\r","").replace(" ","")
    seq = []
    startn =[]
    stopn = []
    for i in range(len(s)):
        if s[i]<='Z' and s[i]>='A':
            seq.append(s[i])
        if s[i]<='9' and s[i]>='0':
            startn.append(s[i])
        if s[i] == '.':
            stopn.append(s[i+2:len(s)])
            break
    return s,''.join(seq),''.join(startn),stopn[0]

def RCSVF(File):
    '''Read a csv file and output its content as a dictionary.Odd columns as keys while even columns as values.'''
    dat= pd.read_csv(File,header = None)
    codon = []
    amino = []
    dictD = {}
    for c in range(np.shape(dat)[1]):
        if c % 2 == 0:
            codon.append(dat.iloc[:,c])
        else:
            amino.append(dat.iloc[:,c])
            for r in range(np.shape(dat)[0]):
                dictD[str(dat.iloc[r,c-1])] = str(dat.iloc[r,c])
    return dictD

def TRANS(SEQ,DICT,s):
    '''Translate a string containing a nucleotide sequence into a string containing the corresponding seqence of amino acids. Allow miss 2 dna pieces at most.'''
    try:
        TRANS = []
        if len(SEQ) % 3 == 1:
            temp = list(SEQ)
            temp.append(SEQ[0:2])
        if len(SEQ) % 3 == 2:
            temp = list(SEQ)
            temp.append(SEQ[0:1])
        if len(SEQ) % 3 == 0:
            temp = list(SEQ)
        if len(temp) % 3 == 0:
            for i in range(0,len(temp),3):
                TRANS.append(str(DICT[''.join(temp[i:i+3])]))
        if len(s) < len(''.join(TRANS))-1:
            print('Miss codon!')
        if len(s) == len(TRANS)-1:
            TRANS = TRANS[:-1]
    except:
        print('condon reading error!')
    return TRANS
        
inputDNA = r"E:\extdata\NM_207618.2.DNA.txt"
inputTRANS =r"E:\extdata\NM_207618.2.trans.txt"
inputTAB =r"E:\extdata\DNA codon Std table.csv"
sDNA = RF(inputDNA)
sTRANS = RF(inputTRANS)
dictTAB = RCSVF(inputTAB)
outputTRANS = TRANS(sDNA[20:937],dictTAB,sTRANS)
# stop word not included(otherwise declude last one using [:-1])
print(''.join(outputTRANS)==sTRANS)

inputFASTA =r"E:\extdata\NM_207618.2.FASTA.txt"
outputFASTA = RFAS(inputFASTA)
sDNA = outputFASTA[1]
startn = int(outputFASTA[2][:])-1
stopn = int(outputFASTA[3][:])-1
sTRANS = RF(inputTRANS)
dictTAB = RCSVF(inputTAB)
outputTRANS = TRANS(sDNA[startn:stopn],dictTAB,sTRANS)