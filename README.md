<h1>Presentation of the project</h1>

TThe main objective of this project was to identify the genomic regions accessible by the DNA during the transcription. In order to do this, we re-analyzed the ATAC-seq results from a publication of 2019 ("STATegra, a comprehensive multi-omics dataset of B-cell differentiation in mouse" David Gomez-Cabrero et al. 2019 https://doi.org/10.1038/s41597-019-0202-7). In this publication, an analysis was performed on a B3 murine cellular line. This cell line from the mouse model the pre-B1 stage. After the nuclear translocation of the transcription factor Ikaros, those cells grow to the pre-BII stage. During this stage, the B cells progenitor are subject to a growth arrest and a differentiation.
This B3 cell line was transduced by a retroviral pathway with a vector coding for a fusion protein, Ikaros-REt2, which can control nuclear levels of Ikaros after exposition to the Tamoxifen drug.
After this treatment, the cultures were collected at t=0h and t=24h.

<h2>Experimental design</h2>

About 50 000 cells were collected for each sample.
There are 3 biological replicates by sample : R1, R2 and R3
Each of the biological replicates was performed for the two cellular stages : 0h and 24h
In total, 6 samples were studied (3 replicates for t=0h and 3 replicates for t=24h)

Each of these samples were sequenced by an Illumina sequencing with Nextera-based sequencing primers. Considering that the sequencing was performed in paired-end, each of the 6 samples has a forward result file and a reverse result file. So, in total, 12 samples were collected for the analysis.

<h3>Raw dataset</h3>

SRR4785152&nbsp;&nbsp;&nbsp;&nbsp;50k_0h_R1_1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.6G&nbsp;&nbsp;&nbsp;&nbsp;50k_0h_R1_2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.6G<br>
SRR4785153&nbsp;&nbsp;&nbsp;&nbsp;50k_0h_R2_1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.5G&nbsp;&nbsp;&nbsp;&nbsp;50k_0h_R2_2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.6G<br>
SRR4785154&nbsp;&nbsp;&nbsp;&nbsp;50k_0h_R3_1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.6G&nbsp;&nbsp;&nbsp;&nbsp;50k_0h_R3_2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.6G<br>

SRR4785341&nbsp;&nbsp;&nbsp;&nbsp;50k_24h_R1_1&nbsp;&nbsp;&nbsp;&nbsp;0.5G&nbsp;&nbsp;&nbsp;&nbsp;50k_24h_R1_2&nbsp;&nbsp;&nbsp;&nbsp;0.5G<br>
SRR4785342&nbsp;&nbsp;&nbsp;&nbsp;50k_24h_R2_1&nbsp;&nbsp;&nbsp;&nbsp;0.5G&nbsp;&nbsp;&nbsp;&nbsp;50k_24h_R2_2&nbsp;&nbsp;&nbsp;&nbsp;0.5G<br>
SRR4785343&nbsp;&nbsp;&nbsp;&nbsp;50k_24h_R3_1&nbsp;&nbsp;&nbsp;&nbsp;0.5G&nbsp;&nbsp;&nbsp;&nbsp;50k_24h_R3_2&nbsp;&nbsp;&nbsp;&nbsp;0.5G<br>

In addition, we must recover the data from previous analysis which are on the HPC cluster. We make this with these commands :

Subsets : scp -r student08@193.49.167.84:/home/users/shared/data/atacseq/data/subset /home/ubuntu/snakemake/data/
Raw data : scp -r student08@193.49.167.84:/home/users/shared/data/atacseq/data/raw /home/ubuntu/snakemake/data/
Bowtie files : scp -r student08@193.49.167.84:/home/users/shared/databanks/Mus_musculus_GRCm39/bowtie2 /home/ubuntu/snakemake/data/reference/
Fasta files : scp -r student08@193.49.167.84:/home/users/shared/databanks/Mus_musculus_GRCm39/fasta/all.fasta /home/ubuntu/snakemake/data/reference/
Trimmomatic tool : scp -r student08@193.49.167.84:/opt/apps/trimmomatic-0.38/trimmomatic-0.38.jar /home/ubuntu/snakemake/
Picard tool : scp -r student08@193.49.167.84:/opt/apps/picard-2.18.25/picard.jar /home/ubuntu/snakemake/
Illumina clips : scp -r student08@193.49.167.84:/home/users/shared/data/atacseq/data/NexteraPE-PE.fa /home/ubuntu/snakemake/data/

<h1>Workflow</h1>
<h2>Quality control of the raw data</h2>

First of all, we need to perform a quality control on the raw data in order to see the global quality of the dataset. In order to do this, we will create an array of 12 elements which contains the sequencing data (forward and rerverse strand). After that, the fastqc function will be used for each of the elements in the array.

<h2>Elimination of the adaptaters</h2>

After the first quality control and before any analysis, we must eliminate the sequencing primers in the sequences in order to not bias the results. So, we use the trimmomatic function on each element of the previous array

The adapters used for the sequencing are Nextera-based primers whose the composition is the following :
>PrefixNX/1
AGATGTGTATAAGAGACAG<br>
>PrefixNX/2
AGATGTGTATAAGAGACAG<br>
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG<br>
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA<br>
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG<br>
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC<br>

<h2>Quality control following the trimming</h2>

After the trimming step, it is necessary to make a second quality control in order to be sure of the quality of the remaining sequences and compare the quality between these sequences and the sequences before the trimming

Once again, we will use the fastqc function but on the results from the trimming.

<h2>Mapping</h2>
Now that our sequences are trimmed, we must synchronize the forward and the reverse file of each sample. We will do this by mapping our sequences on the indexed reference genome (Mus_musculus_GRCm39). The mapping will be done with the Bowtie2 tool.
In order to this, two differnts arrays of 6 elements will be created. In one of them, there will be all the forward files, and in the other there will be the reverse files. The bowtie2 tool, will then take files from each array as an input and map it to the reference genome.

<h2>Elimination of the duplicates</h2>

As a result of the Bowtie2 tool, we get the alignement of the sequences on the reference genome. However, there is some duplicates that can bias the analysis later. That's why we will remove these duplicates thanks to the MarkDuplicates function from the picard module.

<h2>Data exploration</h2>
Now that we have clear results, we can make some statistics on them. These will be done thanks to the Deeptools module. Indeed, this module will permits to know if the samples are correlated (plotCorrelation function).

However, in order to see correlation between the sample, we must combine all the input BAM files into one.

<h2>Identification of DNA access sites</h2>

Thanks to the results from the 'atac_picards.slurm' script we are now able to identify the DNA access sites with the callPeaks function from the MACS2 module. This tool will give results as peaks which will represent the accessible regions.

<h2>Identification of common and unique DNA access sites</h2>

Now that we have the positions of the access regions for the t=0h and the t=24h results, we are now able to see which regions are common and unique between these two stages. This analyze will be performed by the intersect function from the bedtools module. This function will take two files as input and see if there is overlaps between them. 
In order to see the common and unique sites between the two conditions, we need to create two different arrays which contains three elements. One of them will contains all the t=0h files and the other, all the t=24h files.

If an access site is open at t=0h but closed at t=24h, it means that the Tamoxifen drug had an effect on these sites and closed them.
However, if an access site is open at t=0h and is still open at t=24h, it means thant the Tamoxifen drug hadn't an effect on these sites.


<h2>Visualisation with IGV tool</h2>

Now that we know which sites are impacted by the Tamoxifen drug, we can visualize better with the IGV tool.

<h2>Run the program</h2>

You can run all the steps of the workflow with the Snakefile by launching the command :
snakemake --cores all --use-conda --snakefile scripts/Snakefile
<h1>Used tools</h1>

<h2>Fastqc</h2>

fastqc "input_file" -dir "path" -o "path"

  -dir : Path to where the temporary files should be writen<br>
  -o : Path to where the ouput files should be writen<br>

"FastQ Screen: A tool for multi-genome mapping and quality control" Wingett SW, Andrews S. 2018 https://doi.org/10.12688/f1000research.15931.2

<h2>Trimmomatic</h2>

java -jar /opt/apps/trimmomatic-0.38/trimmomatic-0.38.jar PE "input_files" "output_files_names" ILLUMINACLIP: "path" LEADING:"int" TRAILING:"int" SLIDINGWINDOW:"int" MINLEN:"int"

  ILLUMINACLIP : Path to where are the Illumina adaptaters sequence<br>
  LEADING : Cut bases off the start of a read<br>
  TRAILING : Cut bases off the end of a read<br>
  SLIDINGWINDOW : Perform a sliding window trimming<br>
  MINLEN : Specify the minimum length of a read to be kept<br>

"Trimmomatic: a flexible trimmer for Illumina sequence data" Bolger AM et al. 2014 https://doi.org/10.1093%2Fbioinformatics%2Fbtu170

<h2>Bowtie2</h2>

bowtie2  --very-sensitive -p "int" -k "int"  -x "int"  -1 "input_file"  -2 "input_file"

  --very-sensitive : Permits to give better alignment results<br>
  -p : Number of cores to use<br>
  -k : Maximum number of alignment to report per read<br>
  -x : Maximum length of DNA fragment<br>

"Fast gapped-read alignment with Bowtie 2" Langmead, B., Salzberg, S. 2012 https://doi.org/10.1038/nmeth.1923

<h2>Picard - MarkDuplicates</h2>

java -jar /opt/apps/picard-2.18.25/picard.jar MarkDuplicates I= "input_file" O= "path" M= "path" REMOVE_DUPLICATES= "boolean"

  O : Path to where the output files should be writen<br>
  M : Path to where the duplication metrics should be writen<br>
  REMOVE_DUPLICATES : Choose to write the duplicates in the output file with flags or not<br>

“Picard Toolkit.” 2019. Broad Institute, GitHub Repository. https://broadinstitute.github.io/picard/; Broad Institute

<h2>Deeptools</h2>
<h3>MultiBamSummary</h3>

multiBamSummary bins --bamfiles "input_files" -o "path"

  -o : Path to where the output files should be writen

<h3>plotCorrelation</h3>

plotCorrelation -in "input_file" -p "character" --plotNumbers -c "character -o "path" --removeOutliers

  -o : Path to where the output files should be writen<br>
  -p : Choose whether you want to plot a heatmap or a scatterplot<br>
  -c : Choose whether to use the pearson or the spearman coefficient<br>
  --removeOutliers : Will remove outliers from the samples<br>

<h3>bamCoverage</h3>

bamCoverage -b "input_file" -o "path" -of "character"

  -o : Path to where the output files should be writen<br>
  -of : Choose whether the output will be in Bigwig or in bedgraph format<br>

<h3>bamPEFragmentSize</h3>

bamPEFragmentSize -b "input_files" --samplesLabel "character" --table "path"  -o "path"

  -o : Path to where the output files should be writen<br>
  --samplesLabels : Rename the samples by these labels<br>
  --table : Give the metrics in a tabular format<br>

"deepTools: a flexible platform for exploring deep-sequencing data" Ramírez F, Dündar F, Diehl S, Grüning BA, Manke T. 2014 https://doi.org/10.1093%2Fnar%2Fgku365

<h2>MACS2 - CallPeaks</h2>

macs2 callpeak -t "input_files" -f "character" -g "character"  -n "character" --outdir "path"

  -f : Format of the input files<br>
  -n : Name of the output files<br>
  -g : Specify the length of the genome
  --outdir : Path to where the output files shoul be writen<br>

https://github.com/macs3-project/MACS

<h2>Bedtoools - Intersect</h2>

bedtools intersect (-v) -a "input_file" -b "input_file" > "output_path

  -v : Returns entries in -a that have no overlap in -b<br>

"BEDTools: a flexible suite of utilities for comparing genomic features" Quinlan AR, Hall IM. 2010 https://doi.org/10.1093/bioinformatics/btq033
