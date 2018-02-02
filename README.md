# RNA-Seq_analysis


1 Processing of Gene Expression Data
1.1 Sequencing
The RNAseq samples (RNA extraction and sample preparation described previously) for the wild
type H37Rv and the cold sensitive mutant cell lines (2 replicates each) were sequenced at 16 timepoints (∼ 3-hr intervals, 0-55hr) after shifting the temperature from 30◦C (non-permissive) to 37◦C
(permissive). The samples were sequenced on an Illumina 2500 instrument in paired-end mode,
using a read-length of 150+150bp. The mean number of reads per sample was 8.9M(range 4.2-16.5
M). The reads were mapped to the H37Rv genome using the Burrows-Wheeler Aligner (BWA) with
default parameter settings.[1] Reads mapping to each Open Reading Frames (ORF) were totalled
(sense strand only). Counts were truncated to a maximum coverage of 10,000 (reads/nucleotide)
because certain loci were over-represented (e.g. rrs, rnpB, ssr, and Rv3661 had counts 0.5-1M).
1.2 Normalization
Global expression profiles of the cold sensitive mutant samples showed a gradual increase in expression of a few genes that dominate expression at the latter time-points (Figure 1). Consequently, a
compensatory decrease was observed in the expression of other genes, making normalization by the
traditional Reads per Kilo base of transcript per Million (RPKM) mis-representative. Specifically,
3 genes dominated expression in the cold sensitive mutant cells at the latter time-points (30-55hr):
Rv0186A/mymT, Rv0188, and Rv2395A/aprA (Figure 2). To correct for the bias induced by these
outliers, the normalization method implemented in DESeq2 was used, which first normalizes counts
by the geometric mean for each gene across samples, and then scales each dataset to have a common
median (which is less sensitive to outliers).[2] This was applied to all 64 datasets (2 cell lines X 2
replicates X 16 time-points) in parallel. As a result, the expression patterns were well-calibrated
between time-points, with the medians matched (Figure 3). Similarly, a histogram shows that the
distributions of log-fold-changes in the cold sensitive mutant cells relative to the wild type cells were
centred on 0, meaning that on average, genes are equally likely to exhibit an increase or decrease in
expression (i.e. there is no net bias toward increase or decrease over time)(Figure 4). We note that
there is an increase in dispersion in the late time-points.

1.3 Quality Assessment
The two replicates of the cold sensitive mutants and the wild type H37Rv cells show high correlation
between each other indicating good reproducibility (Figure 5 and 6). Similarly, each time-point in
the cos sensitive mutant and the wild type are well correlated with the neighbouring time-points
(see matrix of correlations for replicates among different time-points on spectrum of correlation
coefficients (Figure 7 and 8).The 2-block structure of these diagrams suggest a significant shift in
the expression patterns occurred around 27 hours for the wild type cells and 33 hours for the mutant
cells, probably associated with the first cell division event. Additionally, we observe two sub groups
within the first block for the mutant cells around 18 hrs possibly marking the initiation of cell
division. We note that the time-point 0 for the mutant cells was fairly divergent from all the others,
possibly due to the immediate shock of being release from replication-blockade by temperatureshift; therefore, we dropped time-point 0 from both mutant and wild type samples in the subsequent
analyses.

1.4 Filtering of Genes
To identify a subset of the genes with meaningful expression in the cold sensitive mutant and the
wild type cells, the analysis of variance(ANOVA) test was used along with the false discovery rate
adjustment using the Benhamini-Hochberg procedure. In the cold sensitive mutant, 1105 genes out
of 4018 with ajusted p-value < 0:05 were marked as genes with non significant expression (because
expression patterns for genes with low expression levels with high variance across time-points are
noisy). Figure 9 shows the distribution of the coverage estimates, average expression over all 30
time-points for each gene divided by gene length in nucleotides, for significant and non-significant
genes in the cold sensitive mutant. Here, majority of the non significant genes have 0.25 and lower
coverage estimates, while significant genes have 0.25 and higher coverage estimates. The result shows
that the majority of the non significant genes identified using ANOVA test are the genes with lower
expression values. Similarly, 762 genes were found to be non significant and 3256 genes as significant
in the wild-type cells.
Additionally, a model based clustering was fitted via the Expectation Maximization algorithm. The
optimal number of clusters were identified using the Bayesian Information Criterion to be 21 clusters
among the 2913 significant genes in the cold sensitive mutant.[3] Figure 10 shows the expression profile of the 21 clusters representing decreasing trends (clusters: 1,4,18,20), increasing trends (clusters:
2
9,12,14,21), peaking trends around the middle time-point i.e.∼30hrs (clusters: 8,11,13,15,16,17,19).
2 Operon Analysis of Gene Expression
In order to gain insight into overall data quality and to compare our data with operon models in
literature, we analyzed expression profiles of genes within operons.
We observed that genes within an operon have higher correlation than genes from different operons.
Shell et al. identified 2957 operons in the H37Rv genome based on transcription start site (TSS)
mapping.[4] From this study, we selected 661 operons with two or more genes, and calculated an
average correlation among genes in each operon. On average, the expression pattern of genes in
the same operon were found to have a correlation of 0.41 (cold sensitive mutant) and 0.65 (wild
type H37Rv), while randomly selected genes from different operons showed a correlation of -0.01 to
0.05 (Figure 11). Our results are similar to findings in the literature; for example, Roback et al.
found that the mean correlation of expression for operon pairs and non-operon pairs were 0.60 and
0.16 respectively.[5] Figure 13 and 14 shows an example of ESX-3 operon and rps-rpsl operon with
highly coexpressed genes. Thus, the analysis of gene expression profiles in operons supports that
the quality of our data is good.
Our data provide a global expression profile from synchronized cells that can be utilized to generate
a novel set of operon definitions based on temporal coexpression. We generated an operon model
from our data by sequentially calculating correlation among neighboring genes along the chromosome
(from the same strand with intergenic distance of 100bp or less) and assigning them into operons if
an average correlation among genes is higher than 0.5, which was chosen to maximize consistency
with the operons defined in Shell et al. [1] At a correlation threshold of 0.5, our correlation model
predicted 2971 operons comparable to 2957 operons defined in Shell et al. The model identified 652
operons with two or more genes. Figure 12 shows the distribution of total number of genes in each
operon.
Our model highlights similarities and differences with operon boundaries compared to Shell et al.
[1] For example both our model and that of Shell et al. put eleven genes of ESX-3 secretion
system from Rv0282(eccA3 )-Rv0292(eccE3 ) in a same operon (Figure 13). Conversely, Shell et al.
identified Rv0049-Rv0050(ponA1 ) in a same operon based on TSS but our model separates Rv0049
and Rv0050(ponA1 ) and instead puts Rv0050(ponA1 )-Rv0051 in a same operon (Figure 15). We
3
note that Rv0050(ponA1) and Rv0051 are more likely to be in a same operon as their open reading
frames (ORF) overlap and both of these genes are involved in cell-wall biosynthesis.[6, 7] Generally,
our model captures many operons (esx-3, atp, rps/l, nuo, etc.) with functionally consistent expression
profiles similar to Shell et al.
3 Periodicity Analysis of Gene Expression
A constrained nonlinear curve-fitting method was used to identify genes exhibiting periodic expression in the cold sensitive mutant and the wild type H37Rv. The expression data span a short time
series, which resulted in an expectation of observing periodicity with fractional frequencies such as
1-2 cycles in 55 hours. Common methods for detecting periodicity using a statistical significance
tests such as Fisher’s g test on frequency spectra generated by discrete Fourier transform (DFT)
did not perform satisfactorily on our data because these methods work well only when the Fourier
frequencies of the signal are integral [8, 9, 10].
We fitted the expression values of two replicates to a sinusoidal function using non-linear leastsquares using the Levenberg-Marquardt algorithm.[11] The form of the function used for fitting
was:
y = Asin(wx + φ) + B + Cx
where A is an amplitude, B is a mean offset, w is a frequency, φ is a phase shift, and Cx is a
linear term that allows a gradual increasing or decreasing trend. To obtain reasonable solutions,
the parameters were restricted to the following ranges: amplitude (0∼5), frequency (1∼2), phase
(−π ∼ π), slope (-0.02∼0.02), and mean offset (-2∼2).
The quality of the curve-fitting was determined by the goodness-of-fit measured using the sum of
squares due to error(SSE) (Figure 16). The SSE was calculated as:
SSE =
RP
r=1
NP
i=1
(yir − y ^ir)2
where yi is an observed expression value of a gene at time-point i for replicate r and ^ yir is a predicted
value. A cut-off score of 10.0 was selected empirically based on SSE to identify the genes that fitted
a sinusoidal function with a tolerable variance and an overall periodic trend. The SSE ≤ 10:0 also
corresponds to the correlation coefficient measure ≥ 0:817 between the expression values and the
predicted values in the cold sensitive mutant (Figure 17). The mode of the frequency distribution
was at 1.4, corresponding approximate 30 hrs as the most common cycle time (Figure 18). With
4
an expectation that during 55 hours of data collection, periodic genes should manifest greater than
one and up to two complete cell cycles, we further selected the genes with frequency range 1.3∼2.0
(Figure 18). Genes with frequency closer to lower bounds 1.0∼1.3 were not selected as they displayed
artifacts, sharp curves at the beginning and the end time-points to fit the sinusoidal function. Finally,
genes with amplitude higher than 0.75 were selected as candidate periodic genes.
The analysis using the curve-fitting method identified 380 periodic genes in the cold sensitive mutant
and 161 periodic genes in the wild type H37Rv. This shows that periodic genes are ∼ 2.5x as likely
in the cold sensitive mutant. There were only 26 genes that were found to be periodic in both
mutant and wild type. The distribution of periodic genes for the cold sensitive mutant according
to clusters of orthologous groups(COGs) shows that the genes representing all the COG categories
were present in the periodic group, and the representation for each category is comparable to the
distribution of COG categories among global genes (Figure 19). The periodic genes did not fall
predominantly in any particular functional category. The significantly changed categories in the
periodic group with p-value < 0.05 using Fisher’s Exact Test were transcription (51% enrichment),
post-translational modification (74% enrichment), inorganic ion transport and metabolism (93%
enrichment), and general function (37% reduction). Additionally, a model based clustering was
fitted via the Expectation Maximization algorithm. The optimal number of clusters among the 380
periodic genes in the cold sensitive mutant was determined using the Bayesian Information Criterion
to be six clusters.[3] Figure 20 shows the expression profile of the periodic genes in each of these six
clusters. Note that the clusters vary slightly by frequency (1.39 for cluster 3 to 1.98 for cluster 6),
and also by phase shift (peaking at the different time-points).
The results show that the number of periodic genes in the cold sensitive mutant was significantly
higher than the wild type. Also, the overall variation in the expression profile of genes in the
cold sensitive mutant was substantial compared to the wild type. Similarly, we note that many
genes known to play major roles in replication such as dnaA, dnaK, mprA, mtrA, phoP, and others
were found to be periodic in the cold sensitive mutant (Figure 21). Similarly, figure 22 shows the
expression pattern of all the genes belonging to transcription COG category identified as periodic.
5
4 Temporal peak assignment for gene expression profile
4.1 Smoothing the time series data
We used a Gaussian Process (GP) to integrate and smooth out the raw expression profiles from the
two replicates of the cold sensitive mutants and the wild type H37Rv. A GP model is a Bayesian
model that estimates a probability distribution over functions using Gaussian distributions for likelihood functions.[12] The advantage of a GP is that it is unbiased (does not require assumptions of
form of function), but only assumes that adjacent time-points are weakly coupled together based on
Gaussian distributions.
A Gaussian process is specified by a mean function and a covariance function:
f(x) ∼ GP(m(x); k(x; x0))
A prior mean m(x) = 0 and a covariance function, squared exponential, is given as:
k(x; x0) = σ2 exp(− 1 2
dP
i=1
(xi−x0 i)
l2
i
)
where: l2 = lengthscale; σ2 = variance; d = input dimension
We normalized the expression value e(g; t) (with addition of pesudocounts of 10) of each gene g at
each time-point t by dividing by the mean across all the time-points, and then taking log base e
transformation so that the normalized value e0(g; t) fluctuates with mean of 0. The formula is given
as:
e
0
(g; t) = loge Pe T t(g;t e(g;t ) )
Gaussian estimation of the expression levels for a gene at different time-points, subject to noise, is
given as:
y = f(x) + " where : " ∼ N(µ; σn 2)
The predictive distribution for 15 test time-points (∼ 3-hrs intervals, 3-55 hrs), fx1; x2; :::; x∗g is
specified as:
p(f∗jx∗; x; y) = N(m(x∗); k(x∗))
where:
m(x∗) = k(x∗; x)T (k(x; x) + σ2I)−1y
6
k(x∗) = k(x∗; x∗) − k(x∗; x)T(k(x; x) + σ2I)−1k(x∗; x) + σ2
We utilized the GPy Python package to fit the expression data for each gene independently with zero
mean using the following hyperparameters: variance=1.0, noise variance=0.1, and lenghtscale (range
1 ∼ 50) optimized to the Maximum Likelihood Estimate (MLE) using a grid search method.[13] After
fitting the model, the predicted value (i.e. posterior mean) for each time-point can be extracted.
Figure 23 demonstrates that the fitted values from the GP model generally interpolate between
the observed data at each time-point, but they also present a more smoothed profile by averaging
between adjacent time-points to reduce noise. The error bands in the figure shows uncertainity in
the model.
4.2 Temporal peak assignment
Since most genes in the cold sensitive mutant were not overtly periodic, but did vary in expression
over time, we endeavoured to characterize the time-point at which each gene expression level peaked.
We applied the following algorithmic approach to detect the peak expression time-point for each gene
using the smoothed data obtained from the Gaussian Process analysis (described in earlier section).
The time series T with n observations for each gene with smoothed expression values fx1; x2; :::; xng
at different time-points ft1; t2; :::; tng was defined as:
T = f(t1; x1);(t2; x2); :::;(tn; xn)g
First, to screen out the increasing or decreasing trend at the beginning and the end of a time series,
we excluded the first and last two time-points from the peak assignment. Secondly, to identify
well-spaced major peaks across the time-points, we defined a point xi as a peak if it has greater
magnitude than two nearest neighbours on the both directions.
xi > xi+1; xi+2; xi−1; xi−2 8i = 3;4; :::; n − 2
Furthermore, to filter out the genes with lower fluctuation, the difference between the magnitude of
the highest peak xh and the global minimum gmin was restricted to be greater than 0.5. Additionally,
in case of more than one peak in the time series, all the peaks were constrained to have at least a
half magnitude of the highest peak in the expression profile. Finally, a set of peaks P for a time
series was identified as:
P = f(ti; xi)j(xi > xi+1; xi+2; xi−1; xi−2) ^ (xi − gmin > 0:5) ^ (xi ≥ 0:5 ∗ xh)g 8i = 3;4; :::; n − 2
7
Among the significant genes, the peak assignment identified 1620 genes with a single peak and 71
genes with two peaks in the cold sensitive mutant compared to 903 genes with a single peak and 8
genes with two peaks in the wild type. Similarly, 1211 genes in the mutant and 2344 genes in the wild
type did not have any major peak (Figure 24). The result confirms that the gene expression levels in
the cold sensitive mutant show significantly higher fluctuations than in the wild type highlighting the
effects of synchronous cellular processes in the cold sensitive mutant. Figure 25 shows the clustering
of genes with a single peak in the cold sensitive mutant. Also, we note that among the genes with
a single peak, peak 7 at 27 hours has the highest number of genes, and it includes many genes
associated with cell division and transcription, delineating a probable phase of the cell division in
the cold sensitive mutant (Figure 26).
5 Modeling Gene Expression Data
5.1 Sparse Regulatory Network
The expression level of each gene across 15 time-points can be linearly related with the expression
level of transcription factors.[14] We used a sparse linear modelling approach using the Least Squares
Regression with L1 regularization(Lasso) to model the linear relationship between the cell cycle genes
and their probable regulators i.e. transcription factors.
We model the expression level yit of a gene i at time-point t as a linear combination of the expression
level pjt of transcription factor j at time-point t and an interaction strength aij between gene i and
transcription factor j as
yit =
MP
j=1
pjtaij + it; i = 1; :::; N; t = 1; :::; T
where, N is a total number of genes (non transcription factors), M is a total number of transcription
factors, T is experimental time-points, and it measures error in the model. L1 penalty was applied
on the regression coefficients, an interaction strength, to select the optimal number of regulators for
each gene.
To model the sparse regulatory network, we selected 32 significantly expressed cell cycle genes as
marked by clusters of orthologous groups, and 109 significantly expressed transcription factors identified by Rustad et al.[15] The general housekeeping genes such as sigma factors, chaperonin (groEL),
single stranded binding protein (ssb), antitoxins, and others were filtered from the list of transcription factors .
8
The linear model showed that on average each cell cycle gene was regulated by 13 transcription
factors. Since the number of probable regulators predicted by the model for each gene was fairly
high (Figure 27), we further used the interaction prediction from the linear model to build a network
and identify a basic set of the regulators likely involved in the regulation of cell cycle processes. The
nodes in the network represent either cell cycle genes or transcription factors and the edges represent
interaction between the genes (Figure 28). Similarly, weight of the edges in the network corresponds
to the coefficient of interaction in the linear model.
To obtain the posterior estimate of the network, we used a regulator selection method based on
the greedy algorithm. First, we ranked all the regulators based on the sum of magnitude of coefficients. Second, we used a forward selection approach to select the model with the average coefficient
of determination (R2) of 0.9. The model identified seven transcription factors- Rv3260c/whiB2,
Rv2986c/hupB, Rv2258c, Rv0081, Rv1167c, Rv0880 and Rv0043c as the possible regulators of the
cell cycle genes (Figure 29). The transcription factors Rv3260c/whiB2, Rv2986c/hupB, Rv2258c are
found to positively regulate higher number of the cell cycle genes, while other transcription factors,
Rv0081, Rv1167c, Rv0880 and Rv0043c are found to negatively regulate higher number of the cell
cycle genes with the magnitude of coefficient ≥0.1 in the model (Figure 29b). Also, Figure 30 shows
an example prediction from the network where the cell cylce gene wag31 is positively regulated by
whiB2 and negatively regulated by Rv0043c and Rv1167c transcription f