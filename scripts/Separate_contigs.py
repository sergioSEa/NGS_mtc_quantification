import sys
import gzip
#Input: scripts/Separate_contigs.py $prefix.per-base.bed.gz {output.ncl_depth} {output.mtc_depth} $MTC

Depth_file = sys.argv[1]
Nuclear_depth = sys.argv[2]
Mtc_depth = sys.argv[3]

mtc_contigs = sys.argv[4:]

Nuclear = open(Nuclear_depth, "w")
Mtc = open(Mtc_depth, "w")
with gzip.open(Depth_file, "rb") as DEPTH:
	for position in DEPTH:
		position=position.decode("utf-8") 
		name = position.split("\t")[0]
		if name in mtc_contigs:
			Mtc.write(position)
		else:
			Nuclear.write(position)
Mtc.close()
Nuclear.close()
