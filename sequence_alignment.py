import numpy as np
from helper import Helper

class NeedlemanWuncsh:
	
	def __init__(self, seq1,seq2, gap):
		self.score = 0
		self.seq1 = seq1
		self.seq2 = seq2
		self.gap = gap
		self.blosum62 = Helper.blosum
		
	def matrix_init(self):
		matrix_ = np.zeros( (len(self.seq2)+1, len(self.seq1)+1)) 
		matrix_[0] = [(-i)*self.gap for i in range(len(self.seq1)+1)]
		for i in range(len(self.seq2)+1): 
			matrix_[i][0] = (-i)*self.gap
		return matrix_
	
	def get_value_from_blosum(self,c1,c2):
		v1 = self.blosum62.get((c1,c2))
		v2 = self.blosum62.get((c2,c1))
		if v1 != None:
			return v1
		else:
			return v2
	
	def matrix_nw(self):
		m = self.matrix_init()
		seq1 = self.seq1
		seq2 = self.seq2
		gap = self.gap
		for i in range(1,len(seq2)+1):
			for j in range(1,len(seq1)+1):
				d = m[i - 1][j - 1] + self.get_value_from_blosum(seq1[j-1],seq2[i-1])
				v = m[i - 1][j] - gap
				h = m[i][j - 1] - gap
				m[i][j] = max(d, v,h)
		return m
		
	def get_alignment(self):
		seq1 = self.seq1
		seq2 = self.seq2
		gap = self.gap
	
		m = self.matrix_nw()
		i = len(self.seq1)#kolom
		j = len(self.seq2)#baris
		alignmentSeq1 = ''
		alignmentSeq2 = ''
		alignmentScore = 0
		# while (i >= 0 and j >=0):
		# traceback
		while i+j > 0 :
			score = m[j][i]
			blosum_score = self.get_value_from_blosum(seq1[i-1], seq2[j-1])
			d = m[j-1][i-1] + blosum_score
			h = m[j][i-1] - gap
			v = m[j-1][i] - gap
			
			if score == d:
				alignmentSeq1 = seq1[i-1] + alignmentSeq1
				alignmentSeq2 = seq2[j-1] + alignmentSeq2
				alignmentScore = alignmentScore + self.get_value_from_blosum(seq1[i-1], seq2[j-1])
				i-=1
				j-=1
			elif score == v:
				alignmentSeq1 = '-' + alignmentSeq1
				alignmentSeq2 = seq2[j-1] + alignmentSeq2
				alignmentScore = alignmentScore  -gap
				j-=1
			elif score == h:
				alignmentSeq1 = seq1[i-1] + alignmentSeq1
				alignmentSeq2 = '-' + alignmentSeq2
				alignmentScore = alignmentScore - gap
				i-=1
		a = alignmentSeq1 + '\n' + alignmentSeq2
		return a, alignmentScore
        #print(alignmentSeq1 + '\n' + alignmentSeq2)