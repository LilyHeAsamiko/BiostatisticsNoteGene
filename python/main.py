# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 11:38:01 2024

@author: Qin HeLily 

Run in terminal
"""
import numpy as np
import random
import pandas as pd
from DNAtoolkit import * 
from utilities import colorization

#Creating a random DNAsequence for testing：
DNAStr  = ''.join([random.choice(DNA_neucleotides) for n in range(60)])
inputDNA = r"E:\extdata\NM_207618.2.DNA.txt"
#testStr = read_seqfile(inputDNA)
#DNAStr = ''.join(testStr)
_,testStr,startP,stopP = RFAS(inputDNA)
startPos = int(startP)-1
stopPos = int(stopP)
DNAStr = testStr[startPos:stopPos]
DNA_list = DNAcheck(DNAStr)

print(f'\nSequence: {colorization(DNAStr)}\n')
print(f'[1] + Sequence Length: {len(DNAStr)}\n')
print(colorization(f'[2] + Nucleotide Frequency: {count_DNA(DNAStr)}\n'))

print(f'[3] + DNA/RNA Transcription: {colorization(transcription(DNAStr))}\n')

print(f'[4] + DNA String + Reverse Complement:\n5" {colorization(DNAStr)} 3"')
print(f'   {"".join(["|" for c in range(len(DNAStr))])}')
print(f'3" {colorization(complement(DNAStr))} 5 [complement]"\n')
print(f'5" {colorization(reverse_complement(DNAStr))} 5 [REV. complement]"\n')

print(f'[5] + GC Content: {gc_content(DNAStr)}%\n')
print(f'[6] + GC Content in Subsection k=5: {gc_content_subseq(DNAStr,k = 5)}\n')
#seqtmp = read_seqfile(r'test.fasta')
#print('[6.5] + GC Content in fasta file: {[gc_content(vals) for keys,vals in rosalind_dict(seqtmp).items(])}\n')

#GC content usage one is to detect replication origin()
print(f'[7] + Aminoacids Sequence from DNA: {translate_seq(DNAStr)}\n')
#translated_aminoAcid = translate_seq(DNAStr,0,36)
translated_aminoAcid = translate_seq(DNAStr,startPos,stopPos)
print(f'with frequency: {translated_freq(translated_aminoAcid)}\n')
interested = ['L','R','M']
print(f'[8] + Codon Frequency (L): {codon_usage(DNAStr,interested)}\n')
print(f'[9] + Reading Frame: {ORF(DNAStr)}\n')
print(f'[10] + All Protein Pieces in 6 Open Reading Frames: {all_protein_orfs(DNAStr, 0, len(DNAStr), True)}\n')



