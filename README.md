# Schisto
Scripts used for calculating rare allele sharing and preliminary related clusters

The rare allele sharing script produces output for findSibClusters.pl.

Note: findSibClusters.pl does not actually identify sibling clusters; it does identify clusters of individuals who share proportions of alleles greater than a user-modifiable threshold. Actual realtionships between individuals should be determined by posterior probabilities of different relationship types.

# Example Pipeline
Do variant and sample quality control to create a good-quality vcf.

```perl
perl rareAlleleSharing.pl --vcf myvcf.vcf
```

This will create two files: 
+ The first file contains the proportion of rare alleles shared between all pairwise samples in "list" format. Each row lists two samples in the first two columns and the subsequent columns show the proportion of rare alleles shared between the pair for each generation sampled (default is to sample 500 variants in 30 different generations). Note that allele sharing is directional; the proportion of alleles shared between A and B may be different from the proportion of alleles shared between B and A. This is because the variants are sampled independently between the two comparisons (A to B, then B to A), and because samples A and B may each have different variants that are called and missing. 
+ The second file contains the mean proportion of shared alleles between all samples in matrix format. 

```perl
perl findSiblingClusters.pl --in myvcf.vcf_maxFreq0.1_numVars10_numGens2.rareAlleleSharingPairs > output.txt
```

The above command will create two files:
+ myvcf.vcf_maxFreq0.1_numVars10_numGens2.rareAlleleSharingPairs.sibclst will contain every sample in every putative relative cluster, except for one sample from each cluster. The idea is that this output, which should contain a single column with samples that are putatively closely related to one or more samples in the dataset, could be used to remove these samples from a .vcf.

+ Output.txt will contain  four columns with each line giving information about how a sample is related to a relative cluster in the data. The first column is a relatedness level, which indicates at what level a sample is related to the cluster (see below for a more detailed explanation of this column). The second column is a sample id. The third column is a cluster id, which takes the name of one of the samples in the cluster. Finally, the fourth column is the mean allele sharing of the sample shown in column 2 with members of cluster shown in column 3.

Explanation on “relatedness level”: Lines with a “1” in the first column should indicate “sibling” level relatedness, but this depends on the relatedness thresholds given at runtime (but you might just be using the defaults of 0.45, 0.4, and 0.25). A “2” in the column indicates that the sample is closely related to members of the cluster, but not enough related to be classified as a sibling under default thresholds, and a “3” indicates putative 2nd-degree relatedness to members of the cluster when using default thresholds.
