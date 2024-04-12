# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 11:37:19 2024

@author: QinHeLily
"""
import collections
import pandas as pd
import numpy as  np

DNA_neucleotides = ['A','T','G','C']
DNA_list = {'A':0,'T':0,'G':0,'C':0}
DNA_Complement = {'A':'T','T':'A','G':'C','C':'G'}
inputTAB =r"E:\extdata\DNA codon Std table.csv"

def DNAcheck(seq):
    seq.upper()
    for sq in seq:
        if sq not in DNA_neucleotides:
            return False
    return seq
        
def count_DNA(seq):
    for sq in seq:
        if sq in DNA_neucleotides: 
            DNA_list[sq] += 1
    return DNA_list
    #return dict(collections.Counter(seq))

def transcription(seq):
    '''DNA -> RNA TRanscription. Replacing Thymine with Uracil'''
    return seq.replace('T','U')

def complement(seq):
    '''Swapping adenine with thymine and guanine with cytosine.'''
    return ''.join([DNA_Complement[n] for n in seq])

def reverse_complement(seq):
    '''Swapping adenine with thymine and guanine with cytosine.Reversing newly generated string'''
    return ''.join([DNA_Complement[n] for n in seq][::-1])
    #Pythonic method(a little bit faster)
    #mapping = str.maketrans('ACTG','TAGC')
    #return seq.translate(mapping)[::-1]
    
def gc_content(seq):
    '''GC Content in a DNA/RNA sequence'''
    return round((seq.count('C')+seq.count('G'))/len(seq)*100)

def gc_content_subseq(seq,k = 20):
    '''GC Content in a DNA/RNA subset of a whole sequence'''
    GCpct = []
    for i in range(0,len(seq)-k+1,k):
        subseq = seq[i:i+k]
        GCpct.append(gc_content(subseq))
    return GCpct

def read_seqfile(file):
    with open(file,'r') as f:
        return [lines.strip() for lines in f.readlines()]
    
def rosalind_dict(text):
    rosalind_dict = {}
    for lines in text.readlines():
        if '>' in lines:
            rosalind_dict[lines] = ''
        else:
            rosalind_dict[lines] += lines
#    ['>' in lines?rosalind_dict = '':rosalind_dict[lines] += lines for lines in text.readlines()]
#    [rosalind_dict[lines] = '' if '>' in lines else rosalind_dict[lines] = rosalind_dict[lines]+lines for lines in text.readlines()]
    return rosalind_dict

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

codonDict = RCSVF(inputTAB)

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


def translate_seq(seq, init_pos=0, stop_pos =''):
    '''Translates a DNA sequence into an aminoacid sequence. Maximal length 10000!'''
    try:
        if stop_pos == '':
            stop_pos =min(10000,len(seq))
        assert((stop_pos - init_pos) % 3 == 0)
        assert(stop_pos - init_pos <= 10000)        
        return ''.join([codonDict[seq[pos:pos+3]] for pos in range(init_pos, stop_pos-3+1, 3)])
    except: 
        print('Check sequence translation length(stop_pos - init_pos) requested to be maximal 10000 and devided by 3!')
        
def translated_freq(seq):
    '''output the general frequency of the Translated aminoacid sequence.'''
    tmpDict = dict(collections.Counter(seq))
    freqDict = {}
    for keys,vals in tmpDict.items():
        freqDict[keys] = str(vals/sum(tmpDict.values())*100)+'%' 
    return freqDict

        
def codon_usage(seq,lsts):
    '''count a list of AMINOACID frequencies and its original CODON pieces. '''
    codonLst = []
    aaLst = []
    OutputDict = {}
    for i in range(0,len(seq)-3+1, 3):
        if codonDict[seq[i:i+3]] in lsts:
            codonLst.append(seq[i:i + 3])
            aaLst.append(codonDict[seq[i:i+3]])
    
    freqDict = dict(collections.Counter(aaLst))
    [freqDict.update(e) for e in [{s:0} for s in list(set(lsts)-set(aaLst))]]        
    totalL = sum(freqDict.values())
    for s in lsts:
        if s in aaLst:
            OutputDict[s] = s+':'+str(freqDict[s] / totalL *100)+'% with:'+','.join([codonLst[i] for i,x in enumerate(aaLst) if x == s]) + ' as interested'    
        else:
            OutputDict[s] = s+': 0% as interested'
    return freqDict,OutputDict

def ORF(seq):
    ORF_list = []
    ORF_list.append(translate_seq(seq,0,len(seq)))
    ORF_list.append(translate_seq(seq,1,len(seq)-2))
    ORF_list.append(translate_seq(seq,2,len(seq)-1))
    ORF_list.append(translate_seq(complement(seq),0,len(seq)))
    ORF_list.append(translate_seq(complement(seq),1,len(seq)-2))
    ORF_list.append(translate_seq(complement(seq),2,len(seq)-1))
    return ORF_list

def protein_in_single_orf(seq):
    if '*' in seq and 'M' in seq:
    # protein ending aa == '*'
    # protein starting aa == 'M' 
        return seq[seq.index('M'):seq.index('*')]
    return []

def all_protein_orfs(seq, startPos, stopPos, ordered = True):
    Protein_list = []
    if stopPos > startPos:
        seqTmp = seq[startPos:stopPos]
    else:
        seqTmp = seq
    Protein_list.append(protein_in_single_orf(translate_seq(seqTmp,0,len(seqTmp))))
    Protein_list.append(protein_in_single_orf(translate_seq(seqTmp,1,len(seqTmp)-2)))
    Protein_list.append(protein_in_single_orf(translate_seq(seqTmp,2,len(seqTmp)-1)))
    Protein_list.append(protein_in_single_orf(translate_seq(complement(seqTmp),0,len(seqTmp))))
    Protein_list.append(protein_in_single_orf(translate_seq(complement(seqTmp),1,len(seqTmp)-2)))
    Protein_list.append(protein_in_single_orf(translate_seq(complement(seqTmp),2,len(seqTmp)-1)))
    if Protein_list and ordered == True:
        return sorted(Protein_list,key = len, reverse = True)  
    return Protein_list

