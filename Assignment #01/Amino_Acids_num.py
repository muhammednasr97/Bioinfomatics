import DNA_to_RNA as R

n = len(R.RNA_Seq)
Amino_num, basis = divmod(n, 3)
print('Number of Amino Acids =', Amino_num)
print('Number of basis remaindered =', basis)
