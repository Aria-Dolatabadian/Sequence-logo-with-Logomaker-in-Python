#https://github.com/jbkinney/logomaker/tree/master/logomaker/tutorials

import matplotlib.pyplot as plt
import logomaker as lm
with open('seq.fasta') as f:
    raw_seqs = f.readlines()
print(raw_seqs[10:20])

seqs = [seq.strip() for seq in raw_seqs if ('#' not in seq) and ('>') not in seq]
# preview sequences
print('There are %d sequences, all of length %d'%(len(seqs), len(seqs[0])))
print(seqs[:5])


counts_mat = lm.alignment_to_matrix(seqs)
print(counts_mat.head())
lm.Logo(counts_mat)
plt.show()



#Multiple sequence alignment with gaps

# load ww alignment
with open('ww_sequences.fasta') as f:
    lines = f.readlines()

# preview loaded file
print(''.join(lines[:20]))
# extract ww domain sequences
seqs = [seq.strip().upper() for seq in lines if ('#' not in seq) and ('>') not in seq]

# preview sequences
# Preview sequences
print('There are %d sequences, all of length %d'%(len(seqs), len(seqs[0])))
print(seqs[:10])

# create counts matrix
ww_counts_df = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='.-X')

# preview counts dataframe
ww_counts_df.head()

# show full ww counts
lm.Logo(ww_counts_df)

# filter base on counts
num_seqs = ww_counts_df.sum(axis=1)
pos_to_keep = num_seqs > len(seqs)/2
ww_counts_df = ww_counts_df[pos_to_keep]
ww_counts_df.reset_index(drop=True, inplace=True)

# show cropped ww counts logo
lm.Logo(ww_counts_df)
plt.show()
