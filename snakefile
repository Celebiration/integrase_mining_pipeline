import pandas as pd
import numpy as np
import os
import sys

configfile: "/Yol_Data/NCBI/MAG/Bacteria/refseq_wkdir/configs/config.yaml"
species=int(config["species_taxid"])
output_dir=config["output_dir"]
mapping=pd.read_csv('/Yol_Data/NCBI/MAG/Bacteria/refseq_wkdir/refseq_bacteria_assembly_summary.tsv',sep='\t')

genomes=[i.split('/')[-1] for i in list(mapping.loc[mapping['species_taxid']==species,'ftp_path'])]

if len(genomes)<=1:
    raise ValueError("No more than 1 assemblies found. Exiting.")
else:
    rule all:
        input:
            [output_dir+"/out/"+str(species)+"/"+i+"/per_seq.txt" for i in genomes]

rule gunzip:
    input:
        "{file}.gz"
    output:
        "{file}"
    shell:
        "gunzip {input}"

rule prodigal:
    input:
        "{genome}_genomic.fna"
    output:
        out=output_dir+"/out/"+str(species)+"/{genome}/proteins.gbk",
        protein=output_dir+"/out/"+str(species)+"/{genome}/proteins.faa"
    log:
        output_dir+"/log/"+str(species)+"/{genome}_prodigal.log.txt"
    shell:
        "prodigal -i {input} -o {output.out} -a {output.protein} > {log} 2>&1"

rule find_recombinase:
    input:
        faa=output_dir+"/out/"+str(species)+"/{genome}/proteins.faa",
        fna="{genome}_genomic.fna",
        hmm_file=config["hmm_file"]
    output:
        hmmout=output_dir+"/out/"+str(species)+"/{genome}/hmmout.txt",
        ali=output_dir+"/out/"+str(species)+"/{genome}/ali.txt",
        perseq=output_dir+"/out/"+str(species)+"/{genome}/per_seq.txt",
        #recombinase=output_dir+"/out/"+str(species)+"/{genome}/recombinases.txt"
    threads: 4
    log:
        output_dir+"/log/"+str(species)+"/{genome}_find_recombinase.log.txt"
    shell:
        """
        hmmsearch -o {output.hmmout} -A {output.ali} --tblout {output.perseq} --noali --cut_ga --cpu {threads} {input.hmm_file} {input.faa}
        /Yol_Data/NCBI/MAG/Bacteria/refseq_wkdir/scripts/process_per_seq.sh -f {input.fna} -a {input.faa} -i {output.perseq} -c 0 > {log} 2>&1
        """
