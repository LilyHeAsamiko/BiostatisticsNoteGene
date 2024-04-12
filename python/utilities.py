# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 11:38:01 2024

@author: Qin He Lily

More color sees DEV_DUNGEON
"""
import numpy as np
import random
import pandas as pd
from DNAtoolkit import * 

def colorization(seq):
    bcolors = {
        'A': '\033[92m',
        'C': '\033[94m',
        'G': '\033[93m',
        'T': '\033[91m',
        'U': '\033[91m',
        'reset':'\033[0;0m'
    }

    tmpStr = ''
    
    for n in seq:    
        if n in bcolors:
            tmpStr += bcolors[n] + n
        else:
            tmpStr += bcolors['reset'] + n
    return tmpStr + '\033[0;0m'

def Fibonacci_Loop_Pythonic(number):
    old,new =1,1
    for itr in range(number - 1):
        new, old = old, old + new
    return new

def Fibonacci_Rabbit_Pythonic(months, nbs):
    nParent,nChild =1,1
    for itr in range(months - 1):
        nChild, nParent = nParent, nChild + nbs*nParent
    return nChild+nParent