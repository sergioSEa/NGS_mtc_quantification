#Open some kind of contig file, get numpy/pandas object of contig_name and depth. Perhaps do in R
#Give file for filtering contigs to this script
from Bio import SeqIO
import sys 
from subprocess import Popen, PIPE

fasta = sys.argv[1]
contig = sys.argv[2]
output = sys.argv[3]





def filter_by_length(contig):
	total = 0
	filter_c= []	
	with open(contig, "r") as f:
		for line in f:	
			l = line.split()
			if len(l) != 2:
				continue
			contig = l[0]
			length = int(l[1])
			if length < 2000: #change it was 2000
				filter_c.append(contig)
			else: total += length
		return total, filter_c


total,filter_c = filter_by_length(contig)
filter_c.append("Adam_Phillippy_v6_sli_scf900160275320")
filter_c.append("Adam_Phillippy_v6_sli_scf900160275320")
#CHECK
#with open("../Project/"+organism+"/genome_size.txt", "w") as t:
#        t.write(str(total))

record = list(SeqIO.parse(fasta, "fasta"))
new_list = []	

for r in record:
    if r.id in filter_c: continue
    else: new_list.append(r)

handle = open(output,"w")
SeqIO.write(new_list,handle,"fasta")

  
	
