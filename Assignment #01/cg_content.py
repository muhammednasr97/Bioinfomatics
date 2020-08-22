import math
DNA_Seq = 'ACAGTCGACTAGCTTGCACGTAC'
n = len(DNA_Seq)
c = DNA_Seq.count('C')
g = DNA_Seq.count('G')
cg_percentage = (c + g) * 100 / n
print('The percentage of CG = ', math.floor(cg_percentage), '%')
