from pathlib import Path
from Bio import SeqIO
basedir = workflow.basedir
configfile: "config.yaml"

list_input = config["input_organism"].split(",")
list_nuclear = config["mtc_assembly"].split(",")
list_mtc = config["ncl_assembly"].split(",")
list_forward = config["forward_read"].split(",")
list_reverse = config["reverse_read"].split(",")

dic_targets = {"nuclear":list_nuclear, "mitochondria":list_mtc, "Fread":list_forward, "Rread":list_reverse}

def search(organism, target_file,Inputs= list_input,Targets=dic_targets):
	for n in range(len(Inputs)):
		if organism == list_input[n]:
			out = dic_targets[target_file][n]
			print(out)
			return(out)


			
def trim(reads):
        r = reads.split(".")[0] + "_trimmed." + ".".join(reads.split(".")[1:])    
        return(str(r))



# TARGETS #


rule all:
	input:
		expand('Results/{organism}/metrics',zip,organism=list_input)
	message: "Pipeline complete"


rule contig_size:
	input:  
		ncl = lambda wildcards: search(wildcards.organism, "nuclear"),
		mtc = lambda wildcards: search(wildcards.organism, "mitochondria")
	output:
		ncl = temp("temp/{organism}/ncl_contigs"),
		mtc = temp("temp/{organism}/mtc_contigs")
	resources:
                mem_gb     = 15,
		walltime_h =  1	
	
	run:
		try:
			shell("samtools faidx {input.ncl} ;\n")
		except:
			SeqIO.convert(input.ncl, "fasta", input.ncl, "fasta")
			shell("samtools faidx {input.ncl} ;\n")
		
		try:
			shell(" samtools faidx {input.mtc} ;\n")
		except:
			 SeqIO.convert(input.mtc, "fasta", input.mtc, "fasta")
			 shell(" samtools faidx {input.mtc} ;\n")
		
		shell("cut -f1-2 {input.mtc}.fai > {output.mtc} ;\n" "cut -f1-2 {input.ncl}.fai > {output.ncl} ;\n")


rule filter_nuclear_contigs:
	input:	
		ncl = lambda wildcards: search(wildcards.organism, "nuclear"),
		contigs =  rules.contig_size.output.ncl
	output: 
		temp("temp/{organism}/ncl_2.fa")
	shell:
		"python scripts/filter_contigs.py {input.ncl} {input.contigs} {output} "

rule trimming:
	input:
		R1 = lambda wildcards: search(wildcards.organism, "Fread"),
		R2 = lambda wildcards: search(wildcards.organism, "Rread"),
		adapters = config["adapters"]
	output:
		R1 = temp(trim('temp/{organism}_1.fq.gz')),
		R2 = temp(trim('temp/{organism}_2.fq.gz')),
		stats = "stats/bbduk/{organism}",
		
	params:
		k      = 23,
		min_k  = 11,
		len    = 40,
		mis    = 1,
		qual   = 20,
		poly_g = 10
	threads: 14
	resources:
		mem_gb     = 15,
		walltime_h =  5 # 50Gb WGS takes ~90m
	shell:
		"bbduk.sh in={input.R1} in2={input.R2} out={output.R1} "
		"out2={output.R2} ref={input.adapters} stats={output.stats} "
		"ktrim=r k={params.k} mink={params.min_k} hdist={params.mis} "
		"tpe tbo t={threads} minlen={params.len} qtrim=r "
		"trimq={params.qual} trimpolygright={params.poly_g} "
		"-Xmx5g  \n"

rule mapping:
	input: 
		R1 = rules.trimming.output.R1,
		R2 = rules.trimming.output.R2,
		ncl = rules.filter_nuclear_contigs.output,
		mtc = lambda wildcards: search(wildcards.organism, "mitochondria")
	output:
		reference = 'Results/{organism}/reference.fa',
		bam1 =  temp('Results/{organism}/{organism}_raw.bam'),
		sam = temp('Results/{organism}/{organism}_raw.sam'),
	threads: 36
	resources:
		walltime_h =  5
	shell:
		" cat {input.ncl} {input.mtc} > {output.reference} ; \n"
		" bwa index {output.reference} ; \n"
		"bwa  mem -M -t {threads} {output.reference} {input.R1} {input.R2} > {output.sam} ;\n"
		" samtools view -b -q 30 -@ {threads} {output.sam} |  samtools sort -@ {threads}  -o {output.bam1} ; \n"
rule remove_duplicates:
	input: rules.mapping.output.bam1
	output:
		bam2 =  'Results/{organism}/{organism}.bam',
	resources:  
		walltime_h =  2
	threads: 14
	shell:
		"sambamba markdup -r -t {threads} --tmpdir temp/sambamba/{wildcards.organism} {input} {output.bam2} ;\n"
rule depth:
	input: 
		bam = rules.remove_duplicates.output.bam2,
		contig_mtc = rules.contig_size.output.mtc
	output:
		mtc_depth = 'Results/{organism}/depth_mtc',
		ncl_depth = 'Results/{organism}/depth_ncl',
	resources:
		mem_gb     = 120,
		walltime_h =  10
	threads: 3
	shell: 
		"samtools index {input.bam} ;\n"
		"prefix=Results/{wildcards.organism}/depth ; \n"
		"MT=$(awk '{{print $1}}' {input.contig_mtc}) ; \n"
		"mosdepth -t {threads} $prefix {input.bam} ; \n"
		"python scripts/Separate_contigs.py $prefix.per-base.bed.gz {output.ncl_depth} {output.mtc_depth} $MT ; \n"
rule average_depth:
	input: 
		depth_file = 'Results/{organism}/depth_{sample}'
	output:
		computation_result = 'Results/{organism}/Estimation_{sample}'
	resources:
		mem_gb = 100
	shell:
		"Rscript scripts/calculate_statistics.r {input.depth_file} {output.computation_result} ;\n"

rule metabolic_estimation:
	input:
		ncl_depth = 'Results/{organism}/Estimation_ncl',
		mtc_depth = 'Results/{organism}/Estimation_mtc',
	output:
		output_file = 'Results/{organism}/metrics',
		conc = temp('Results/{organism}/depth')
	run:
		shell(" cat {input.mtc_depth} {input.ncl_depth} > {output.conc}  ;\n ")
		coverage = []
		ratio = 0
		#Median Depth \t Coverage
		with open(str(output.conc),"r") as f:
			for lines in f:
				if lines[0] == "M": continue
				lines = lines.rstrip()
				lines = lines.split(",")
				coverage.append(lines[1])
				if ratio == 0: ratio = float(lines[0]) 
				else: ratio = ratio/float(lines[0])
		s = [wildcards['organism'],str(ratio)]
		s.extend(coverage)
		with open(output.output_file,"w") as out:
			out.write("#Organism\tRatio of Median Depth\tCoverage_mtc\tCoverage_ncl\n")
			out.write("\t".join(s)+"\n")
