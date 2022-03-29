# coding=utf-8
from Bio import SeqIO
from baseFun import *
'''
计算每个位置上ATCG的使用率
'''



file = 'E.coli/IMemProt.cds'
# file ='B.subtulis/B.subtulis.cds'
seqList = readFasta(file, 'list')
print len(seqList), 'genes'

f = open('E.coli/nuclRate.txt', 'w')
# f = open('B.subtulis/nuclRate.txt', 'w')
f.write('first\t\t\t\tsecond\t\t\t\tthird\n')
f.write('A\tT\tC\tG\tA\tT\tC\tG\tA\tT\tC\tG\tlength\n')
for seq in seqList:
	nuclDict = getCondonNuclotideSeq(seq[:-3])
	for key in [1,2,3]:		# 每个位置的核苷酸
		nuclSeq = nuclDict[key]
		N = float(len(nuclSeq))
		A, T, C, G = nuclSeq.count('A'), nuclSeq.count('T'), nuclSeq.count('C'), nuclSeq.count('G')
		f.write('%.6f\t%.6f\t%.6f\t%.6f\t' % (A/N, T/N, C/N, G/N))
	GC = (seq.count("G")+seq.count("C"))/float(len(seq))
	f.write('%d\n' % len(seq))
f.close()




