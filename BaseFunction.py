# coding=utf-8
from Bio import SeqIO

def readFasta(file, flag):
	'''
	read a fasta file and return the records as a dict or list
	by id as key and seq as value
	'''
	fastaDict = {}
	seqList = []
	terminatorList = ['TAG', 'TGA', 'TAA']
	hande = SeqIO.parse(file, 'fasta')
	for record in hande:
		ID = record.id
		Seq = str(record.seq).upper()
		if len(Seq)%3!=0 or Seq[-3:] not in terminatorList:		# 进行筛选
			print ID
			continue
		if flag=='dict':
			fastaDict[ID] = Seq.upper()
		elif flag=='list':
			seqList.append(Seq.upper())

	if flag=='dict': return fastaDict
	elif flag=='list': return seqList

def median(dataList):
	'''
	caculate the meidan of the dataList
	the result will be a float or int
	'''
	data = sorted(dataList)
	size = len(data)
	if size%2==0:
		median = (data[size/2] + data[size/2-1]) / 2.0
	elif size%2==1:
		median = data[(size-1)/2]
	return median

def classSeqByLength(seqList):
	'''
	devide the seqList as longer and shorter group
	the mid length seq will be devided to shorter group
	return these two group as a dict
	where keys are 'longer' and 'shorter', and seq as value
	'''
	classedSeq = {}
	sortedSeq = sorted(seqList, key=lambda x:len(x))
	size = len(seqList)
	if size%2==0:
		midIndex = (size/2)-1
	elif size%2==1:
		midIndex = (size-1)/2
	classedSeq['part1'] = sortedSeq[midIndex+1:] 	# long seq
	classedSeq['part2'] = sortedSeq[0:midIndex+1]	# short seq
	return classedSeq

def devideSeq(seq):
	'''
	devide codons as front and latter parts
	the mid codon will be devided to latter part
	'''
	length = len(seq)
	midIndex = (length/3/2)*3
	front = seq[0:midIndex]
	latter = seq[midIndex:]
	return front, latter

def getCondonNuclotideSeq(seq):
	'''
	get nuclotide in each codon position
	return the nuclotides as three seq
	'''
	nuclSeq = {1:'', 2:'', 3:''}
	for i in range(0, len(seq)-2, 3):
		nuclSeq[1] += seq[i]
		nuclSeq[2] += seq[i+1]
		nuclSeq[3] += seq[i+2]
	return nuclSeq


def countNuclotide(seq):
	'''
	count the nuclotide number in seq
	return a list
	ordered by A T C G
	'''
	nuclotideNumber = {}
	nuclotideNumber[A] = seq.count('A')
	nuclotideNumber[T] = seq.count('A')
	nuclotideNumber[C] = seq.count('A')
	nuclotideNumber[G] = seq.count('A')
	return nuclotideNumber


def getFounctionalGeneLocation(pttFile):
	'''
	find funtional genes in ptt file
	return a locationList
	ptt_location: 190..255
	'''
	f = open(pttFile)
	f.readline()
	f.readline()
	f.readline()
	locationList = []
	for line in f.readlines():
		if 'hypothetical' in line or 'putative' in line or 'unknown' in line\
		or 'uncharacterized' in line:
			continue
		locationList.append(line.split('\t')[0])
	return locationList


def getSeqInffn(locationList, ffnFile):
	'''
	mach gene sequence in ffn by location
	ffn_location: 190-255
	'''
	seqList = []
	for record in SeqIO.parse(ffnFile, 'fasta'):
		ID = record.id
		seq = record.seq
		temp = ID.split(' ')[0]
		temp = temp.split(':')[-1]
		if 'c' in temp:		# 使 ffn_location 与 ptt_location 格式一致
			ffn_location = temp.split('-')[1] + '..' + temp.split('-')[0].replace('c', '')
		else:
			ffn_location = temp.replace('-', '..')
		if ffn_location in locationList:
			seqList.append(seq)
	return seqList

