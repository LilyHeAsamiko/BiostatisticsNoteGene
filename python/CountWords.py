import pandas as pd
import numpy as np

def RF(File, decluders):
    '''Read a txt file and output its content as a string, replacing decluders as spaces.'''
    f =open(File, "r")
    s= f.read()
    for d in decluders:
        s=s.replace(d," ")
    return s

def CountWords(P, spliter = ' ', Capital = True):
    '''Read a piece of words and count the words. With spliter as space, difference between capital and small defaulltedly '''
    if Capital == False:
        P = P.lower()
    if spliter == ' ':
        P=P.split()
    else:       
        for d in spliter:
            P=P.split(d)    
    wordsCount = {}
    for word in P:
        if word in wordsCount:
            wordsCount[word] += 1
        else:
            wordsCount[word] = 1
    return wordsCount

def CountWordsF(P, spliter = ' ', Capital = True):
    '''Read a piece of words and count the words. With spliter as space, difference between capital and small defaulltedly '''
    if Capital == False:
        P = P.lower()
    if spliter == ' ':
        P=P.split()
    else:       
        for d in spliter:
            P=P.split(d)    
    wordsCountF = Counter(P)
    return wordsCountF
        
inputf = r"C:\Users\Admin\Documents\CountWordstest.txt"
Ptest = RF(inputf,[',','.','"',"'",':',';','?','!','\n','\t','\r'])
WCD = CountWords(Ptest)
print(WCD)
from collections import Counter
WCDF = CountWordsF(Ptest)
print(WCD==WCDF)
WCD = CountWords(Ptest).find('I')