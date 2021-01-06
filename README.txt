--------------------------------------------------
-------		INSTALLATION 		----------
--------------------------------------------------

The IMAGE pipeline is distributed as ready-to-use code, which means there is no 
installation per se. It is only necessary to inflate the compressed archive, 
which you also had to do to look into this file. So congratulations on installing IMAGE.

Add the IMAGE folder to $PATH for easy execution. If you recieve errors or IMAGE doesn't execute,
please ensure that IMAGE.pl and all files in the utils folder are executable on your system.

--------------------------------------------------
------		DEPENDENCIES		----------
--------------------------------------------------

Running the IMAGE pipeline depends on having the statistical software R
installed and executable with the several packages (glmnet, methods, cluster, 
doParallel, foreach, Matrix, data.table, GenomicRanges) installed.

Please see http://www.r-project.org/ for instructions on installation of R.

Please see http://www.bioconductor.org/ for instructions on package installation
within R

Additionally, IMAGE uses HOMER (homer.salk.edu). 
Please make sure that HOMER is installed and fully executable.

Alternatively, you can use conda to quickly install all necessary dependencies 
from the environment.txt file

---------------------------------------------------
------		INPUT			-----------
---------------------------------------------------

IMAGE requires three input files: 

One contains the full genome fasta sequence. This file can be downloaded from various sources, e.g. the UCSC genome browser.

One contains gene expression data. This file is most easily generated using iRNA-seq.
Otherwise, a raw tag count matrix (with header line) with the following columns is needed:
RefSeq	Symbol	Chr	Start	End	Strand	Blank_Column	Gene_Length_For_RPKM_Calculation	Raw_Count_A	Raw_Count_B	...

The other file contains enhancer data. This file should have the following columns (without a header):
Chr	Start	End	PeakID	Normalized_Count_A	Normalized_Count_B	...

Run IMAGE.pl for an explanation of command-line options needed in IMAGE.

---------------------------------------------------
------		OUTPUT			-----------
---------------------------------------------------

IMAGE produces an .R file as output. This file contains several data.frames
and lists with the results:

GeneRPKM: Contains RPKM normalized gene expression data and differential expression analysis.
Enhancers: Contains the supplied enhancer information 
Target_Genes: List of data.frames of each motif with target genes (indicated by Target column; 1 = Target)
Target_Sites: List of data.frames of each motif with target sites (indicated by Target column; 1 = Target)
Enhancer_Activity: Contains motif activity of all motifs for each sample using only enhancers (IMAGE step 1)
Gene_Activity: Contains motif activity of all motifs for each sample (Full IMAGE model) 
Result: Contains the processed results. The CausalTF row indicates hits (1 = high confidence, 2 = medium confidence).
	WeightedPvalue: IMAGE uses this value to select candidates. It is a combined score using changes in motif acticity, change in gene expression and maximal expression level.
	Evidence: Indicates the evidence type behind the motif. 
		Direct: Direct assay of the transcription factor (eg. ChIP-seq, SELEX). 
		Indirect: Determined using an indirect assay. 
		Inferred: Inferred by IMAGE. 
		ZifRC: Inferred using ZifRC.
	Pearsons: Contains the Pearsons correlation coefficient between the expression pattern of the transcription factor, and the motif activity of the associated motif. 
	MinimalActivityPvalue:	The lowest p-value of testing motif activities; all-vs-all conditions.

---------------------------------------------------
------		TESTING			-----------
---------------------------------------------------

Example files of gene expression data and enhancer activity
is available in the examples folder. 

IMAGE.pl -region examples/Enhancers.txt -expression examples/GeneExpression.txt -fasta PATH_TO_MM9_GENOME_FASTA -RNADesign 1 1 2 2 3 3 4 4 -EnhancerDesign 1 1 2 2 3 3 4 4 -p 12 -n Test

Should produce a file called Test.R in your current directory demonstrating the pipeline, and the following output in your terminal:

#### Welcome to IMAGE ####

Your region file contains 10000 regions and 4 conditions across 8 files
Your exon file contains 5000 genes and 4 conditions across 8 files
Your fasta file is /data/Genomes/mouse/mm9/mm9.fa
You are using 12 processors.

Starting the analysis
        Setting up for parallizing motif search
        Preparing input file for motif searching
        Scanning for motifs
        Running analysis in R
                1. Reading data into R
                2. Analyzing gene expression
                3. Converting motif hits to matrix - Takes a while
                4. Calculating motif activity - Full model stage
                5. Performing ridge regression
                        Sample 1 completed
                        Sample 2 completed
                        Sample 3 completed
                        Sample 4 completed
                        Sample 5 completed
                        Sample 6 completed
                        Sample 7 completed
                        Sample 8 completed
                6. Calculating motif activity - Reduced model stage
                7. Creating weight matrices
                8. Calculating factor activity - Full model stage
                        Sample 1 completed
                        Sample 2 completed
                        Sample 3 completed
                        Sample 4 completed
                        Sample 5 completed
                        Sample 6 completed
                        Sample 7 completed
                        Sample 8 completed
                9. Calculating factor activity - Reduced model stage
                10. Processing results
                11. Outputting results


IMAGE completed in 362 seconds

 
