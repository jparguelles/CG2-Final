# CG2-Final

## ABOUT

This pipeline is meant to handle the assembly, annotation and downstream analyses of spider silk and venom sequences. Ideally, this pipeline will allow for the assembly of problematic genome data, as well as targeted annotation of silk and venom genes using curated BLAST databases. These sequences are then separated from the remaining genome data for downstream analysis, which can occur parallel to the annotation of remaining sequences through BLAST analyses vs the nr database. Silk sequences are then passed through a python script which identifies and quantifies amino acid motif data, outputting files which can be used in graphical analysis. Venom proteins are moved through additional annotation steps to further identify various structural motifs important for inferring protein activity. 

Beginning from raw long read data (in .fastq format), the pipeline will offer some initial QC via NanoPlot (long-read QC), followed by assembly via Flye (which is purported to better handle assembly of highly repetitive data, as expected in spider genomes). The resulting assembled genomes will then be assessed for completeness via BUSCO (using the arachnid database), before moving on to annotation. The initial annotation step is threefold, with the first step handled via an existing, in-house pipeline (Tr-TAP) which is specifically tailored to spider silk sequences. For this project, we will focus solely on the initial Tr-TAP step, which is focused on BLAST analysis. (Additional steps of the Tr-TAP pipeline can also be incorporated, if desired, and are focused on cleaning up and filtering the resulting data.) The Tr-TAP pipeline is then adapted to venom peptides, through the incorporation of curated BLAST databases which are substituted for the existing silk databases. Sequences not annotated as either silk or venom, are separated out and then further investigated through a BLAST search using the NCBI nr database. Finally, the silk and venom proteins are taken through a few additional analyses. Silk genes are run through a python script (written by me and to be cleaned up to be less ugly eventually), which identifies and quantifies various amino acid motifs associated with mechanical performance of fibers. Venom proteins are taken through additional annotation steps for structural annotation, to identify known structural motifs which can be used to infer their potential function. 

Some additional steps can be done, dependent upon available data. Specifically, should short-read DNA reads be available, a short-read error correction can be performed on the assembled genome prior to gene annotation. Similarly, in the case gland-specific transcriptomic data is available, these genes can be cross-referenced with recovered silk and venom genes from the genome assembly. While not strictly necessary for overall venom/silk gene analysis, this can be used to filter out those genes which are being actively transcribed venom and silk genes… which are perhaps more relevant when assessing questions of venom/silk protein functional properties. Further, annotation via Tr-TAP can be adapted to suit alternative systems, using a personally curated database in place of the databases presented in this document. Similarly, the motif identification step allows for the user to designate amino acid motifs other than those built-in, so can be adapted to seek out other motifs of interest in sequence data.

### Required Programs
NanoPlot (https://github.com/wdecoster/NanoPlot)
FiltLong (https://github.com/rrwick/Filtlong)
Flye (https://github.com/fenderglass/Flye)
BWA (https://github.com/lh3/bwa)
PolyPolish (https://github.com/rrwick/Polypolish)
BUSCO (https://busco.ezlab.org/)
BUSCO Databases: https://busco-data.ezlab.org/v4/data/lineages/
BLAST (https://www.ncbi.nlm.nih.gov/books/NBK279690/)
Tr-TAP (https://github.com/thclarke/TrTAP; private)
RSEM, used in Tr-TAP pipeline: https://github.com/deweylab/RSEM 
SeqTK (https://github.com/lh3/seqtk)
InterProScan (https://github.com/ebi-pf-team/interproscan)


# Step 1: Concatenation and initial QC of long reads

This step operates under the assumption that the long read data being used is in .fastq format. If raw reads exist in multiple files, they can be concatenated into a single .fastq.gz file before moving on to assembly (zipping is not strictly necessary, but will reduce file size):

```
cd WORKINGDIRECTORY
cat *fastq > SPECIES_ALL_READS.fastq
gzip SPECIES_ALL_READS.fastq
```

Once all files to be assembled are concatenated, quality control can be performed using NanoPlot, a program specifically tailored to long-read data:

```
NanoPlot –threads NUMBERTHREADS –verbose –store –tsv_stats \
--N50 –loglength –outdir OUT_DIR –fastq \
SPECIES_ALL_READS.fastq.gz
```

While most flags (--vebose, --store, --tsv_stats –N50 –loglength) in the above code relate to specific QC measures to be visualized in the .html output, others will need specific input from the user. These flags are denoted via ALLCAPS input.

The main result of this analysis will be a .html output which can be viewed to gain insights into the overall quality of read data. This report will provide more data than is necesarrily needed, so below are summarized some key takeaways:

**In Summary Statistics:** 
n50 relates to the overall length of reads in your input data. For this analysis, a higher number (and thus longer read lengths) are ideal. This data is further summarized in histograms lower on the page. It is worth looking over the bell curve (weighted histogram after log transformation) in order to decide whether it makes sense to filter out smaller reads based on overall length distributions. Smaller reads, particularly those below 1kb are unlikely to be very informative for the assembly, so can potentially be filtered out via FiltLong:

```
filtlong –min_length LENGTH SPECIES_ALL_READS.fastq.gz > SPECIES_1kb.fastq
```

An additional step can also be used to discard the lowest quality reads (10% in this example) from those remaining:

```
filtlong –keep_percent 90 SPECIES_1kb.fastq
gzip SPECIES_1kb.fastq
```

## Step 2: Genome Assembly

Flye is chosen as the assembler for this pipeline, as it is supposedly better at handling very repetitive sequence data. As spider silk genes are extremely repetitive, this is desirable:

```
 cd WORKDIR

flye --nano-raw SPECIES_1k.fastq.gz -g ESTIMATEDGENOMESIZE -o OUTDIR -t NUMTHREADS --asm-coverage 50
```

**Some caveats with the above code are:** 

--nano-raw is used if the input data is ONT read data. If you are assembling pacbio data, you would use the –pacbio-raw flag instead.

-g is not necessary as of Flye v2.9, so can be excluded UNLESS you are also using the --asm-coverage flag, which specifies an additional filtration step only necessary when dealing with particularly large genomes. If this does not apply, you may exclude both flags.

--resume can also be used should your run quit mid-analysis for whatever reason (such as running out of storage or memory), and will allow for the run to pick up where it left off

### OPTIONAL STEP: Short Read Polishing

Should you also have short read DNA data available, you can use it here to polish the Flye output (which will be named assembly.fasta, if not renamed). This step utilizes both bwa-mem2 and the python script ‘polypolish_insert_filter.py’:

```
cd WORKDIR

bwa-mem2 index assembly.fasta
bwa-mem2 mem -t NUMTHREADS -a assembly.fasta ../SHORT_READS1.gz > ALIGN1.sam
bwa-mem2 mem -t NUMTHREADS -a assembly.fasta ../SHOT_READS2.gz > ALIGN2.sam
```

```
polypolish_insert_filt.py \
--in1 ALIGN1.sam --in2 ALIGN2.sam \
--out1 filt1.sam --out2 filt2.sam
```

```
polypolish assembly.fasta filt1.sam filt2.sam > ASSEMBLY_POLISHED.fasta
```

## STEP 3: Assembly Quality Control

As before the assembly, it is important to assess the quality of the resulting Flye assembly. As spiders often lack a reference genome, we will focus on looking completeness via BUSCO (using the arachnid database):

```
cd WORKDIR

busco -i ASSEMBLY_POLISHED.fasta -o ASSEMBLY_BUSCO -l LOCALPATHTOARACHNIDDB -c NUMTHREADS -m geno -f --out_path PATHTOOUTDIR –offline
```

While most flags require specific input from the user, others require minor explanation: 

‘-m geno’ denotes that the input is genomic data

‘--offline' forces the run to circumvent attempts to pull databases directly from the BUSCO servers… and while this may sound appealing, for whatever reason attempting to run “online” seems to stall out and kill runs, so do not forgo this flag!

Off the output, we’re most interested in the .txt filed titled ‘short_summary…’ which will look something like this:

	***** Results: *****

	C:90.2%[S:74.9%,D:15.3%],F:6.9%,M:2.9%,n:1013	   
	914	Complete BUSCOs (C)			   
	759	Complete and single-copy BUSCOs (S)	   
	155	Complete and duplicated BUSCOs (D)	   
	70	Fragmented BUSCOs (F)			   
	29	Missing BUSCOs (M)			   
	1013	Total BUSCO groups searched		   

## STEP 4: Annotation

As spider silk genes can be particularly problematic, our lab has developed an in-house pipeline for dealing with the sequences called ‘Tr-TAP’. (THIS PIPELINE IS NOT PUBLIC, SO I CAN’T SHARE GITHUB PAGE, but am including relevant submission scripts, etc)

While the pipeline itself also includes steps for cleaning up and filtering read data, considering things like chimeric sequences and transposable elements, this document will focus on the BLAST submission step for initial annotation. 

This step utilizes curated BLAST databases. In this case, they are databases containing spider silk genes, but this can easily be substituted for a curated database of some other gene-set using the following scripts: 

```
makeblastdb -in FASTAFILE -dbtype NUCL/PROT 
```

(NUCL for DNA/RNA sequences, PROT for protein sequences)

After the database is generated, Tr-TAP can be used to rapidly parse the assembled genome for genes contained in the curated database:

```
cd WORKDIR

perl run_blast_submit.pl -t ASSEMBLY_POLISHED.fasta -b PATHTODATABASE -o OUTDIR -r -f RSEMDIR -i OUTHEADER -c Fly -k -q -p NUMTHREADS
```

Most flags are user-specified, though some may need additional explanation:

-f specifies the directory created for generation of fastq files during a paired-end RSEM run

-I denotes the output header (such as species name)

-c denotes a specific database (in this case, the Fly database) to run a chimera test on, can be excluded

-k skips BLAST vs nr database, if nr BLAST is desired this flag can be excluded

-q auto-generates PBS submissions during the run, so can be excluded if not using a PBS submission system

This will automatically run BLAST analyses vs all databases contained within the database directory, as well as create RSEM databases for additional analysis.

This step can then be repeated on any number of databases, to generate various subsets of BLAST annotations from the assembled genome. In this case, two separate runs would be done, one using silk sequence databases and an additional run using venom sequences databases. 

If access to Tr-TAP is unavailable, you can run independent BLAST analyses for each database generated using the following script:

```
blastx -query POLISHED_ASSEMBLY.fasta -db PATHTODATABASEDIRECTORY -out OUTDIR -num_threads NUMTHREADS -evalue 1e-5 -max_target_seqs 1 -outfmt 6
```


**NOTE**: ‘blastx’ is used here, which automatically translates DNA/RNA sequences to protein sequences. This is not recommended for large query datasets, or for BLAST analysis vs large databases (such as nr), as it attempts 6 BLAST analyses per query sequence. As an alternative, you can substitute ‘blastn’ to BLAST nucleotide sequences or, if query sequences have already been translated into protein sequences, ‘blastp’.


Additionally, some flags can be altered:

-evalue denotes a MAXIMUM evalue for cutoff, this can be changed to the user’s preference

-max_target_seqs denotes the number of BLAST results to report in the output. As this pipeline uses *very* curated databases, the first output is often sufficient. However, if BLASTing against larger databases (such as NCBI’s nr database) it is recommended to request additional output, as the first hit is not necessarily the best.

-outfmt denotes the desire output format and can be specficied to the user’s desired output according to: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

Ultimately, this will result in two subsets of data: genes annotated as spider silk sequences and those annotated as venoms.

### OPTIONAL STEP: Annotation of Remaining Genes

A subset of genes in the assembled genome, minus those annotated as silks and venom can be further annotated via BLAST analysis versus the nr database. This can be run in parallel with other analyses, in the background, if desired. This step, however, requires the removal of previously annotated genes from the assembly, so as not to waste time re-annotating. This can be done via SeqTK’s ‘subseq’ process. 

**NOTE**: *.fasta is used if you have multiple outputs to combine into a single gene list, otherwise, use the name of the file you’d like to generate the gene list from

First, use SeqTK to generate a fasta file consisting SOLELY of your silk/venom gene output. This will be useful for additional analyses detailed later in this document, so keep it handy. 

This begins by generating a list of all annotated genes via their IDs. As the BLAST output can be easily opened as a table (via Excel, for example), it’s possible to simply copy and paste these IDs into a text file. Alternatively, this can also be done via awk:

```
awk ‘{print $1}’ *.fasta > SEQTKLIST.lst
```

Once this list has been generated, you can use SeqTK to grab these genes and put them into a separate .fa, venom example below:

```
seqtk subseq POLISHED_ASSEMBLY.fasta SEQTKLIST.lst > SEQOUTVENOM.fa
```

For here, you can leverage this new list to REMOVE these sequences from the overall assembly. First, filter the venom/silk lists out of the assembly file:

```
grep -f SEQTKLIST.lst POLISHED_ASSEMBLY.fasta -v > ASSEMBLY_REMAINING.lst
```

Then use seqtk as above to assembly a subset of sequences, consisting only of the remaining sequences:

```
seqtk subseq POLISHED_ASSEMBLY.fasta ASSEMBLY_REMAINING.lst > ASSEMBLY_NOVENOM.fa
```

From this point on, the silk and venom genes will undergo separate analyses, as the data hoping to be gleaned from them is different. 

**NOTE**: If transcriptomic data is available, you may also want to consider focusing on just the genes which are being transcribed in venom/silk gland specific transcriptomic data. This can be accomplished by BLAST-ing the transcriptomic data against the assembled genome (POLISHED_ASSEMBLY.fasta converted into a database) and identifying overlap. This would follow the same procedure described above for generating custom BLAST databases.

## STEP 5: Amino Acid Motif Identification and Quantification in Full-Length Silk Genes

For genes identified as silk in previous steps, and containing both N- and C- termini, we can use a python script to generated output consisting of identification, counts and positions of motifs commonly associated with silk fiber mechanical properties. This (goofy and non-optimized) script also allows for in input of desired amino acid motifs, if you are for some reason not studying spider silk genes. It is easiest run in jupyter notebook :

(code provided separate doc, it’s very ugly I swear I’ll optimize it eventually)

## STEP 6: Structural Annotation of Venom Genes

As much of the interest in venom genes comes from better understanding their modes of action, it is important to further annotate these genes by looking at their overall structures. For this document, we’ll focus on functional annotations via InterProScan. 

**NOTE**: Tr-TAP automatically incorporates some functional annotation as part of the existing pipeline, so if using Tr-TAP, this steo in not necessary but can be useful for comparing results. Similarly, as Tr-TAP automatically generates protein sequences as part of its output, this step operates under the assumption that you have translated your gene data.

In the event your data has been translated through some process other than Tr-TAP, it is important to makes sure there are no hidden characters:

```
sed ‘s/*//g’ SEQTKOUTVENOM.fa > SEQVENCLEAN.fasta
```

These files should already be small enough to pass through InterProScan, but in the event you’ve used larger databases for annotation, it is useful to split up large fasta files and parallelize the InterProScan run (uses .pl script from class, so not openly available):

```
cd WORKDIR

splitFASTA.pl -i SEQVENCLEAN.fasta -o OUTDIR -s NUMSEQS
```

**NOTE**: -s denotes the number of sequences per output file. As an alternative, you can use the ‘-f’ flag to instead denote total number of files to be created.

From here, you can run files through InterProScan (this also uses in-class script):

```
cd WORKDIR

interproscan.sh --input SEQVENCLEAN.fasta --seqtype p --formats TSV,XML,GFF3 --goterms --iprlookup --cpu 8 --pathways --disable-precalc --disable-residue-annot --excl-applications MobiDBLite,Phobius-1.01 --output-dir OUTDIR --tempdir TEMP
```




