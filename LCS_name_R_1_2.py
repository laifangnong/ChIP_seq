import numpy as np
def longest_common_subsequence(name1,name2):
	len1=len(name1)
	len2=len(name2)
	matrix=np.zeros((len1,len2))
	for i in np.arange(len1):
		for j in np.arange(len2):
			if name1[i]==name2[j]:
				if i!=0 and j!=0:
					matrix[i][j]=int(matrix[i-1][j-1])+1
				else:
					matrix[i][j]=int(1)
			elif i!=0 and j!=0:
				matrix[i][j]=max(int(matrix[i-1][j]),int(matrix[i][j-1]))
	#print(np.max(matrix))
	return(len1-int(np.max(matrix)))
#longest_common_subsequence('123456','123446')
#longest_common_subsequence('abc','abcd')
#longest_common_subsequence('x1_na_R1.fastq','x2_na_R2.fastq')
