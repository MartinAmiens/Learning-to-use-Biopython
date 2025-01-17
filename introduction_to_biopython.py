from Bio.Seq import Seq
from Bio import SeqIO
my_seq = Seq("AGTACACTGGT")
my_seq
print(my_seq)
my_seq.complement()
my_seq.reverse_complement()

for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    
for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))