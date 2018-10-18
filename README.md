# Comparison of Clustering Methods for Single-Cell RNA Sequencing Data

##### By Ji-Eun Park, Yitong Wang, Davide Risso

## Table of Contents
[1. Introduction](#Introduction)

[2. Results](#Results)

[3. Discussion](#Discussion)

[4. Methods](#Methods)

[5. References](#Reference)

### **Introduction** <a name="Introduction"></a>

Single cell RNA sequencing(scRNA-seq) is widely used to profile the transcriptome of individual cells. From each cell, mRNA is isolated and reverse-transcribed to cDNA for high-throughput sequencing. The number of reads mapped to each gene is then used to quantify its expression in each cell. This allows biological resolution that cannot be matched by bulk RNA sequencing, at the cost of increased technical noise and data complexity. enables detailed profiling of different cell populations and can be used to reveal lineage relationships or discover new cell types.

Sequencing is the technique in measuring RNA. In the earlier stage of RNA sequencing, the RNA was measurable only in a bulk level where the features of each single cell were unable to be distinguished. This was called Bulk RNA sequencing. The technology to measure RNA in a single cell level was rather newly found. In 2009, there were only one high-throughput sequencing of RNA from single cell measured. As the technology improved, now we are able to gain datasets with a size greater than one million cells and with better quality.

In reality, cell types are still unknown many RNA cells. The advantage of clustering single cell RNA, is that by clustering the dataset or cells into different groups, we could predict the unknown cell types based on the known cell types in the same group, or even discover new, unique cell types.

There are many proposed methods in clustering methods on single cell RNA sequencing data. Many methods have been tested, developed, and validated using simulated datasets. However, until now there are no consensus in which method suits best for this particular data type. There are several reasons to this. First, current simulations are often poorly documented, the similarity to real data is not demonstrated, or no reproducible code is available. Second, in reality we do not know the true number of clusters which makes difficult to apply certain methods, such as K-means or PAM, and also makes it difficult to conclude whether a clustering method is clustering the dataset well.

Thus, our goal is to compare the results of different clustering methods using simulations of real datasets as well as real datasets which have the actual cell types included so that we can measure how well each method works. Two simulated datasets and one real dataset were used in this project and the quality of the methods were calculated by the average rand index.

### **Results** <a name="Results"></a>

In the project, for each simulated or real dataset, principal component analysis(PCA) was done and used to perform the two clustering methods, K-means and Partitioning Around Medoids(PAM). The methods were compared using the adjusted Rand Index (aRI) calculated with the clustering results and the simulated or the real group. The method with a higher adjusted Rand Index is considered to be a better method. Details regarding to the simulation process, PCA, and adjusted Rand Index can be found in **Methods**.

There were three simulated datasets using the first actual dataset from the study on analyzing the single-cell RNA-seq data from differentiating HBC stem cells of the olfactory epithelium (Simulation 1). For each simulation, K-means and PAM was performed and was compared to the groups assigned in the simulation. Clustering works better with K-means algorithm in general given a sequence of different number of principal component and the number of clusters  (figure 1 left panel). It can also be seen that the sequential methods produces more stable aRIs compared to the non-sequential or traditional methods (figure 1 right panel). However, the highest average random index was found in the non-sequential methods when the number of clusters were set to the true number of clusters, which in this case was three groups.

![Picture1](https://ws3.sinaimg.cn/large/006tNbRwly1fwc5tt6eh7j30r307f0ut.jpg)

> **Figure1** Simulation results of the first dataset (p63-HBC-diff). Left panel of the first graph shows the average rand index by the number of principal components(pc) for non-sequential K-means and PAM. Right panel of the first graph shows the average rand index by the number of principal components(pc) for sequential K-means and PAM. Red lines represent K-means and the navy lines represent PAM. The second graph shows the average rand index by the number of centroids or medoids (K) for both non-sequential and sequential K-means and PAM. The left panel of the second graph using pc=5 and the right panel using pc=10. The solid lines represent non-sequential methods and the dashed lines represent sequential methods.



The second dataset was also simulated into three different simulation datasets and then the same process was done (Simulation 2). Similar to the results of the simulation of the first dataset, clustering works better K-means algorithm works consistently better in general given a sequence of different number of principal component (figure 2 left panel). However, it also shows that when the correct number of principal components and the number of true clusters were selected, the average Rand Index calculated by the number of clusters is higher for PAM (figure 2 right panel). This is inconsistent to the result of simulation 1, where K-means worked better than PAM. In addition, the sequential method seems to be noisy.

![Picture2](https://ws4.sinaimg.cn/large/006tNbRwly1fwc5txjwz7j30r307fab8.jpg)

> **Figure 2** Simulation results of the second dataset (human colorectal tumors). Left panel of the first graph shows the average rand index by the number of principal components(pc) for non-sequential K-means and PAM. Right panel of the first graph shows the average rand index by the number of principal components(pc) for sequential K-means and PAM. Red lines represent K-means and the navy lines represent PAM. The second graph shows the average rand index by the number of centroids or medoids (K) for both non-sequential and sequential K-means and PAM. The left panel of the second graph using pc=5 and the right panel using pc=20. The solid lines represent non-sequential methods and the dashed lines represent sequential methods.

Unlike the first and second dataset, the third dataset was not simulated and the original dataset. The dataset was processed with the same steps except for simulation. Unlike the simulation results, there is no clear superior method between K-means and PAM given a sequence of different number of principal component (figure 3 left panel). There seems to be no difference in the accuracy of the two methods when
the average Rand Index was calculated by different number of clusters (figure 3 right panel). Hence, unlike simulated datasets, it seems that in reality, the best clustering method is not dominated by neight K-means nor PAM.

![Picture3](https://ws1.sinaimg.cn/large/006tNbRwly1fwc5u0w11sj30r307f75e.jpg)

> **Figure 3** Results of the third dataset (primary glioblastoma) without simulation. Left panel of the first graph shows the average rand index by the number of principal components(pc) for non-sequential K-means and PAM. Right panel of the first graph shows the average rand index by the number of principal components(pc) for sequential K-means and PAM. Red lines represent K-means and the navy lines represent PAM. The second graph shows the average rand index by the number of centroids or medoids (K) for both non-sequential and sequential K-means and PAM. The left panel of the second graph using pc=5 and the right panel using pc=10. The solid lines represent non-sequential methods and the dashed lines represent sequential methods.

![Picture4](https://ws3.sinaimg.cn/large/006tNbRwly1fwc5u3p45dj30r307fmyr.jpg)

> **Figure 4** Plots of first and second principal component of dataset 3. Panel A is colored by the true clusters or cell types given in the real dataset. Panel B is colored by the K-means clustering results. Panel C is colored by the PAM clustering results. Colors are randomly assigned, where the same color does not mean that they are actually the same group or cell type. PAM and K-means work similarly, but does not perfectly match the real clusters result.

### **Discussions**<a name="Discussions"></a>

In general, for the simulated datasets, K-means algorithm works better than Partitioning Around Medoids(PAM) algorithm when true number of principal components and the number of clusters were correctly assigned. However, it is even seen that for some cases, when the number of clusters was correct, PAM works better than K-means. Moreover, when applied to a real data, neither of the algorithms significantly dominates as a superior method no matter which number of principal component and cluster are chosen.

From this difference in the results of simulations and the real data, it can be seen that it is difficult to simulate realistic data. This might be caused by the way the simulation model is designed in the *splatter* package. The model is designed to simulate from the mean of the clusters. To be more specific, the model specifies the mean and simulates the data around the mean, which makes it better for K-means which use the centroids than PAM which means medoids. Hence, it would be difficult to figure out a method that is superior when applied to real datasets even though they are found to be working better on simulations. 

Despite the limit of simulation, it is shown from the simulation that the non-sequential or traditional clustering algorithms gives the highest average random Index. However, this only true when the number of clusters (k) is correctly assigned. The highest average random index of sequential methods were never higher than the non-sequential methods. The advantage of sequential algorithms, however, is that the average RIs are robust to the choice of the number of clusters which coincides with the purpose of sequential algorithms. The downside of these methods is that they take approximately twice as much time compared to traditional algorithms and requires a higher computational power. 

As a conclusion of this project, there were no perfect method in our datasets. There are still many possible proposed clustering methods can be applied such as subsampling, ensemble clustering, etc. However, in order to find out the best clustering method, the primary objective would be to find a simulating algorithm that enables to generate a realistic simulation. The reason is in real life, it is more often the case that the real cluster, group, or cell types are unknown and the number of datasets with the real cluster is minimal. In order to find a consensus on the best algorithm, it is necessary that we simulate realistic data to apply and verify the optimal method.  In the future, a number of combinations of parameters can be applied to the simulation algorithm from the *splatter* package. 

### **Methods**<a name="Methods"></a>

#### Simulated datasets:

*splatter* package was used to simulate our two datasets, p63-HBC-diff dataset(2017) and Li dataset(2017). The source of those two datasets can be found in **Data Availability & Description**. “The core of the Splat model is a gamma-Poisson distribution used to generate a gene by cell matrix of counts. Mean expression levels for each gene are simulated from a gamma distribution and the Biological Coefficient of Variation is used to enforce a mean-variance trend before counts are simulated from a Poisson distribution. Splat also allows you to simulate expression outlier genes (genes with mean expression outside the gamma distribution) and dropout (random knock out of counts based on mean expression). Each cell is given an expected library size (simulated from a log-normal distribution) that makes it easier to match to a given dataset.” (Zappia, 2018).  The package *splatter* can be installed by 

```
source(https://bioconductor.org/biocLite.R); 
biocLite("splatter")
```

The simulation involves the following steps:

1. Set up simulation object
2. Simulate library sizes
3. Simulate gene means
4. Simulate groups/paths
5. Simulate BCV adjusted cell means
6. Simulate true  counts
7. Simulate dropout
8. Create final dataset

The final output is a “[SingleCellExperiment](http://127.0.0.1:36953/help/library/SingleCellExperiment/html/SingleCellExperiment.html)” object that contains the simulated counts but also the values for various intermediate steps. These are stored in the [colData](http://127.0.0.1:36953/help/library/splatter/help/colData) (for cell specific information), [rowData](http://127.0.0.1:36953/help/library/splatter/help/rowData) (for gene specific information) or [assays](http://127.0.0.1:36953/help/library/splatter/help/assays) (for gene by cell matrices) slots.

The parameters which were varied in the simulations are `group.prob`,` method`, `verbose`, `de.facLoc`, `de.facScale`. Details on the parameter choices follow.

`group.prob`: The probabilities that cells come from particular groups. We set probability of cells coming from particular groups to be 0.3, 0.3, 0.4

`method`: which simulation method to use. Options are "single" which produces a single population, "groups" which produces distinct groups (eg. cell types) or "paths" which selects cells from continuous trajectories (eg. differentiation processes). Here we use method equals “groups”

`verbose`: set the verbose argument to FALSE to stop Splatter printing progress messages   

`de.facLoc`: Location (meanlog) parameter for the differential expression factor log-normal distribution.

`de.facScale`: Scale (sdlog) parameter for the differential expression factor log-normal distribution.

Taking p63-HBC-diff data as an example, firstly, we omit missing values and omit the rows who has all 0 in the row and use the cleaned dataset to generate parameters by function `splatEstimate()`, then we use the following code to simulate 3 different data. 

```
lapply(1:3, function(i) splatSimulate(params,group.prob = c(0.3,0.3,0.4), method = "groups",verbose = FALSE,de.facLoc=.25, de.facScale=1,seed=i))
```

#### **Real datasets:** 

We use the real Patel dataset from human tumor tissue. Patel et al. (2014) examined the genome sequence of single cells isolated from brain glioblastomas. This dataset has 5948 features, 430 samples and researcher cluster this real dataset into 9 levels of cell type. You could find the source of Patel et al. in **Data Availability & Description**.                 



#### **Normalization methods:**  

If you look at the column sums of the counts, you can see that each sample had a different amount of reads which could be aligned to the different genes. So if we’re going to do some analysis, we definitely need to take care of the fact that the different samples had different sequencing depth. “It’s not a great idea to use the number of reads which map to genes as a normalizing factor. There have been a number of papers showing that, this number is susceptible to being altered by a single gene. You might have one or a dozen genes which have very high gene expression, and those genes play the largest role in determining this sum. There are a number of more robust estimators. Probably using any of those robust estimators is better than using this sum.” (Sthda.com, 2018) We normalize data in order to scale raw count values to account for the “uninteresting” factors. By doing so, the expression levels are more comparable between and/or within samples. Also, We use logcount to reduce the influence of outliers, one way to make data distribution more like Gaussian data to make PCA works best. As normalization is essential, PCA is applied to normalized account. 

As for our dataset, firstly, we normalized the simulated data and then logcount the normalized date and next we sort out rows that are not all 0 logcounts.



#### **Dimensionality reduction methods:** 

We use PCA as our dimensionality reduction methods and it is applied to p63-HBC-diff, Li and Patel dataset.

**PCA**: According to《An Introduction to Statistical Learning》, "principal components analysis (PCA) is a technique for reducing the dimension of a (n × p) data matrix 'X'. The first principal component direction of the data is that along which the observations vary the most.” and “Principal component analysis (PCA) refers to the process by which principal components are computed, and the subsequent use of these components in understanding the data. PCA is an unsupervised approach, since it involves only a set of features X_1 , X_2, . . . , X_p, and no associated response 'Y'. Apart from producing derived variables for use in supervised learning problems, PCA also serves as a tool for data visualization" (James et al., n.d.)

We used the function`runPCA` from the R package *scater* for the simulation study and, for computational efficiency for the real data sets. We are not sure which number of PCs are appropriate, so we specify five different values of PCs, 2,5,10,20,50.  We used the following parameters: `ncomponents = 50, method = "irlba", ntop = 1000.`



#### **Clustering methods:**

We use four different clustering methods to assess the dataset, K-means, PAM, sequential K-means and sequential PAM.             

##### **1. K-Means:** 

The goal of this algorithm is to find groups in the data, with the number of groups represented by the variable 'K'. The algorithm works iteratively to assign each data point to one of 'K' groups based on the features that are provided. Data points are clustered based on feature similarity(Trevino, 2018). Given a set of observations X = {X1, X2,..., Xn}, where each observation is a d-dimensional real vector, k-means clustering aims to partition the n observations into k (<= n) sets S = {S1, S2, …, Sk} so as to minimize the within-cluster sum of squares (WCSS) (i.e. variance). Formally, the objective is to find:         

<p align="center">
  <img src="https://ws3.sinaimg.cn/large/006tNbRwly1fwc5xh1ayhj30cf03c74h.jpg">
</p>


##### **2. PAM:** 

PAM stands for “partition around medoids”. The algorithm is intended to find a sequence of objects called medoids that are centrally located in clusters. Objects that are tentatively defined as medoids are placed into a set 'S' of selected objects. If 'O' is the set of objects that the set 'U = O − S' is the set of unselected objects. The goal of the algorithm is to minimize the average dissimilarity of objects to their closest selected object. Equivalently, we can minimize the sum of the dissimilarities between object and their closest selected object (Cs.umb.edu, 2018).​                                       

We used the following parameters for Non-sequential K-means and PAM in `clusterMany` function: 

```
clusterMany(ks=2:20, alphas=0.1, betas=0.8, minSizes=1, clusterFunction= c("kmeans","pam"), sequential = FALSE, subsample=FALSE, reduceMethod="PCA", nReducedDims =  c(2,5,10,20,50), verbose=TRUE)
```

All the other parameters were left at their default values. Since we do not know the number of clusters, so we cluster with K-means and PAM with k from 2 to 20 and take only number of components of interest.

##### **3. Sequential K-means and Sequential PAM:** 

Another way to modify the k-means/PAM procedure is to update the means one example at a time, rather than all at once. Clustering approaches are implemented in the `clusterMany` function of the `clusterExperiment` package. Difference in the code is we add the parameter `sequential=TRUE`

```
clusterMany(alphas=0.1,betas=0.8,minSizes=1, clusterFunction=c("kmeans","pam"), sequential=TRUE, subsample=FALSE, reduceMethod = "PCA", nReducedDims = c(2,5,10,20,50), verbose=TRUE)
```



#### **Evaluation criteria:**                                                              

In order to compare clustering results against external criteria, a measure of agreement is needed. Since we assume that each gene is assigned to only one class in the external criterion and to only one cluster, measures of agreement between two partitions can be used.(Faculty.washington.edu, 2018). We use Adjusted Rand Index instead of Rand Index. The reasons are explained below.

As for Rand Index, Given a [set](https://en.wikipedia.org/wiki/Set_(mathematics)) of 'n' [elements](https://en.wikipedia.org/wiki/Element_(mathematics)), S = {o1, o2,.., on} and two[ partitions](https://en.wikipedia.org/wiki/Partition_of_a_set) of 'S' to compare, X = {X1, X2,.., Xr\}, a partition of 'S' into 'r' subsets, and Y = {Y1,Y2,..,Ys}, a partition of 'S' into 's' subsets, define the following:

- a: Number of pairs of elements that are same subset in X and in same subset in Y
- b: Number of pairs of elements that are different subsets in X and in same subset in Y
- c: Number of pairs of elements that are same subset in X and in different subsets in Y
- d: Number of pairs of elements that are different subsets in X and in different subsets in Y



The Rand index, R, is:

<p align="center">
  <img src="https://ws4.sinaimg.cn/large/006tNbRwly1fwc5xtf0tmj306802ga9w.jpg">
</p>

Intuitively, `a+b` can be considered as the number of agreements between X and Y and `c+d` as
the number of disagreements between X and Y.

Since the denominator is the total number of pairs, the Rand index represents the frequency of occurrence of agreements over the total pairs, or the probability that X and Y will agree on a randomly chosen pair.

A problem with the Rand index is that the expected value of the Rand index of two random partitions does not take a constant value (say zero). The adjusted Rand index proposed by (Hubert and Arabie, 1985) assumes the generalized hypergeometric distribution as the model of randomness. 

Adjusted Rand Index is the corrected-for-chance version of the Rand index.

<p align="center">
  <img src="https://ws4.sinaimg.cn/large/006tNbRwly1fwc5svkke0j30jw05lq3j.jpg">
</p>

where n_ij, a_i, b_i are values from the contingency table.

We implemented in the `comparing.Partitions` function of the *clusterSim* package. 



#### **Data Availability & Description:**

- p63-HBC-diff (<https://github.com/rufletch/p63-HBC-diff> )
  - Fletcher RB*, Das D*, Gadye L, Street KN, Baudhuin A, Wagner A, Cole MB, Flores Q, Choi YG, Yosef N, Purdom E, Dudoit S, Risso D, Ngai J. Deconstructing Olfactory Stem Cell Trajectories at Single Cell Resolution. Cell Stem Cell (2017)
  - 47083 features, 849 samples

- Li dataset (<https://hemberg-lab.github.io/scRNA.seq.datasets/human/tissues/> )
  - Li, H. et al. Reference component analysis of single-cell transcriptomes elucidates cellular heterogeneity in human colorectal tumors. Nat. Genet. (2017)
  - 55186 features, 561 samples

- Patel dataset (<https://hemberg-lab.github.io/scRNA.seq.datasets/human/tissues/>)
  - Patel, A. P. et al. Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. Science 344, 1396–1401 (2014)
  - 5948 features, 430 samples



### **References**<a name="References"></a>

1. Datta S, Nettleton D, (2014) *Statistical Analysis of Next Generation Sequencing Data,* Springer International Publishing
2. Risso D, 2014, Cluster Analysis, lecture notes,  Division of Biostatistics PB HLTH C240D/STAT C245D, University of California, Berkeley, delivered 02/12/2014
3. Zappia L, Phipson B, Oshlack A (2017). “Splatter: simulation of single-cell RNA sequencing data.”*Genome Biology*. doi:[ ](http://doi.org/10.1186/s13059-017-1305-0)[10.1186/s13059-017-1305-0](http://doi.org/10.1186/s13059-017-1305-0),[ ](http://dx.doi.org/10.1186/s13059-017-1305-0)<http://dx.doi.org/10.1186/s13059-017-1305-0>.
4. Purdom E, Risso D (2018). *clusterExperiment: Compare Clusterings for Single-Cell Sequencing*. R package version 2.0.2.
5. Faculty.washington.edu. (2018)[ ](http://faculty.washington.edu/kayee/pca/supp.pdf)<http://faculty.washington.edu/kayee/pca/supp.pdf>  [Accessed 23 Jul. 2018]
6. Sthda.com. (2018). *RNA sequencing data analysis - Counting, normalization and differential expression - Easy Guides - Wiki - STHDA*. [online] Available at: http://www.sthda.com/english/wiki/rna-sequencing-data-analysis-counting-normalization-and-differential-expression#why-rnaseq-data-should-be-normalized [Accessed 1 Aug. 2018].
7. James, G., Witten, D., Hastie, T. and Tibshirani, R. (n.d.). *An introduction to statistical learning*. Pp.230,231,375.
8. Trevino, A. (2018). *Introduction to K-means Clustering*. [online] Datascience.com. Available at: https://www.datascience.com/blog/k-means-clustering [Accessed 2 Aug. 2018].     
9. Cs.umb.edu. (2018). [online] Available at: https://www.cs.umb.edu/cs738/pam1.pdf [Accessed 2 Aug. 2018].                 



