# epicv2refine

<b>This function deals with replicate/duplicate probes on the EPICv2 array in two ways.</b> For each group of replicates:
1. Probes mapped to "chromosome 0" in the Illumina manifest, and probes labelled as inferior (due to reduced sensitivity and precision) by <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-024-10027-5">Peters et al., BMC Genomics 2024</a>, are removed. Often this leaves only one 'superior' replicate probe in the group.
2. If there is more than one probe left in the group after removing inferior probes*, mean and variance are calculated across the remaining probes. For each sample, beta values of remaining replicate probes are replaced by mean of beta values. E.g.</p>
Before:</p>

| Probe Name | Sample 1 | Sample 2 |
| -------- | -------- | -------- |
| cg00002033_TC11 | 0.6533183 | 0.4650523 |
| cg00002033_TC12 | 0.6440894 | 0.4638761 |

</p>
After:  </br>

| Probe Name | Sample 1 | Sample 2 |
| -------- | -------- | -------- |
| cg00002033_TC1_12 | 0.6487038 | 0.4644642 |

</p>

<b>High-variance probes</b></br>
By default, if probes in the replicate group have unusually high variance in beta values compared to other replicate groups (outliers defined as variance > Q3 + 1.5*IQR, same as in a box plot) then all probes in the replicate group will be removed. Alternatively, user can set remove_high_var = FALSE, to keep and average these high-variance replicates as shown above.
  
*<i>Reasons for this include that one probe may have the best sensitivity while another has the best precision, group mean may be superior to any one probe, or probes may have no inferior/superior label due to insufficient evidence.</i>

