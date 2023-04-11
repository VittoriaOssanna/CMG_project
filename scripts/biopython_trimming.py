from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from collections import Counter

inp_aln = AlignIO.read('core_gene_alignment.aln', 'fasta')
thr = 0.875 # 1 out of 8
cols_to_remove = set()
sub_aln = []

for i in range(len(inp_aln[0])):
  nrows = len(inp_aln) - inp_aln[:, i].count('-')

  for count in Counter(inp_aln[:, i].replace('-', '')).values():
    if (count / nrows) >= thr:
      cols_to_remove.add(i)
      break
   
cols_to_remove_set = set(cols_to_remove )

for aln in inp_aln:
  seq = ''.join([c for i, c in enumerate(aln.seq) if i not in cols_to_remove])
  
  if seq:
    sub_aln.append(SeqRecord(Seq(seq), id=aln.id, description=''))

  with open('core_gene_alignment.trim.aln', 'w') as f:
AlignIO.write(MultipleSeqAlignment(sub_aln), f, 'fasta')
