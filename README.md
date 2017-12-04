# Integrated analysis of motif activity and gene expression changes of transcription factors
**Jesper G. S. Madsen<sup>1,4</sup>, Alexander Rauch<sup>1,4</sup>, Elvira Laila Van Hauwaert<sup>1</sup>, Søren Fisker Schmidt<sup>1,2</sup>, Marc Winnfeld<sup>3</sup>, Susanne Mandrup<sup>1,5</sup>**

<sup>1</sup> Department of Biochemistry and Molecular Biology, University of Southern Denmark, Odense, Denmark.<br>
<sup>2</sup> Present address: Institute for Diabetes and Cancer, Helmholtz Center Munich, German Research Center for Environmental Health, Neuherberg, Germany.<br>
<sup>3</sup> Research and Development, Beiersdorf AG, Hamburg, Germany.<br>
<sup>4</sup> These authors contributed equally.<br>
<sup>5</sup> Corresponding author (s.mandrup@bmb.sdu.dk)<br>

### Abstract
The ability to predict transcription factors based on sequence information in regulatory elements is a key step in systems-level investigation of transcriptional regulation. Here, we have developed a novel tool, IMAGE, for precise prediction of causal transcription factors based on transcriptome profiling and genome-wide maps of enhancer activity. High precision is obtained by combining a near-complete database of position weight matrices (PWMs), generated by compiling public databases and systematic prediction of PWMs for uncharacterized transcription factors, with a state-of-the-art method for PWM scoring and a novel machine learning strategy, based on both enhancers and promoters, to predict the contribution of motifs to transcriptional activity. We applied IMAGE to published data obtained during 3T3-L1 adipocyte differentiation and showed that IMAGE predicts causal transcriptional regulators of this process with higher confidence than existing methods. Furthermore, we generated genome-wide maps of enhancer activity and transcripts during human mesenchymal stem cell commitment and adipocyte differentiation, and used IMAGE to identify positive and negative transcriptional regulators of this process. Collectively, our results demonstrate that IMAGE is a powerful and precise method for prediction of regulators of gene expression. 

### Scripts for reproduction of manuscript figures
All figures (except those not prepared using R) in the manuscript can be reproduced from the R scripts within this repository.<br>
To reproduce the figures locally, please clone this repository and download the data files [here](http://bioinformatik.sdu.dk/solexa/webshare/IMAGE/IMAGE_Data.tar.gz).<br>
Each script in markdown format (.md) can be opened in RStudio and the code chunk within backticks can be executed.<br>

Main figures
-------------
[Figure 1](Links/Figure1.md)<br>
[Figure 2](Links/Figure2.md)<br>
[Figure 3](Links/Figure3.md)<br>
[Figure 4](Links/Figure4.md)<br>
[Figure 5](Links/Figure5.md)<br>

Supplementary figures
-------------
[Figure S1](Links/FigureS1.md)<br>
[Figure S2](Links/FigureS2.md)<br>

### Links to relavant sites
[IMAGE download](http://bioinformatik.sdu.dk/solexa/webshare/IMAGE/IMAGE_v1.1.tar.gz)<br>
Full-length manuscript<br>
[NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104537)<br>
[Mandrup group website](http://sdu.dk/mandrupgroup)<br>
