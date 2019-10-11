# NGS_mtc_quantification
Snakemake pipeline for mitochondrial quantification from NGS reads.  


### Prerequisites

* python3.  
 pathlib.  
 Biopython.  
 Snakemake.   

* R  
 dplyr  
 readr.  
 ggplot2  
* bwa
* samtools
* bbduk
* sambamba
* mosdepth


### Installing

Once all the requirements are in order,  clone the repository

```
git clone https://github.com/sergioSEa/NGS_mtc_quantification.git
```

### Running

Prepare the config file (config.yaml) with the following paths:  

* mtc_assembly: path to the location of the assembly of mitochondrion in FASTA format.  
* ncl_assembly: path to the location of the assembly of the nuclear genome in FASTA format.  
* forward_read: Path to the location of forward reads (.fq), which can be both gzipped (.gz) or not. 
* reverse_read: Path to the location of forward reads (.fq), which can be both gzipped (.gz) or not.  
* input_organism: This is the ID use in all output files
* adapters: Path to adapter sequences in FASTA format, used for trimming.  

Multiple files can be specified at the same time so they will be run in parallel on a cluster. Example:
```
* mtc_assembly: Fasta1,Fasta2  
* ncl_assembly: Fasta_nuclear1,Fasta_nuclear2  
* forward_read: forward1,forward2 
* reverse_read: reverse1,reverse2  
* input_organism: ID1,ID2
```
Run by activating the smk file with Snakemake

```
snakemake -s MetaQuantNGS.smk -np
```
This will run a dry run, i.e, the pipeline will be displayed but not run. If no errors are reported, the pipeline might be run removing the dry-run tag

```
snakemake -s MetaQuantNGS.smk
```

Please refer to Snakemake documentation for automatic submission of jobs to diffrent clsuter queuing systems.

## Authors

* **Sergio Andreu-Sanchez** - *Initial work* - [sergioSEa](https://github.com/sergioSEa)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


