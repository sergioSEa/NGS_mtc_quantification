from Bio import SeqIO
import sys

Name = sys.argv[1]#"scaffold_4794"

assembly = sys.argv[2] #"Data/ncl_assemblies/B10K-DU-009-16.genomic.fa"


record = list(SeqIO.parse(assembly, "fasta"))
new_list = []

for r in record:
    if Name == r.id: continue
    else: new_list.append(r)

handle = open(assembly,"w")
SeqIO.write(new_list,handle,"fasta")


