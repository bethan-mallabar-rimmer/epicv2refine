# epicv2refine

<b>This function deals with replicate/duplicate probes on the EPICv2 array in two ways.</b> For each group of replicates:
1. Probes labelled as inferior (due to reduced sensitivity and precision) by <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-024-10027-5">Peters et al., BMC Genomics 2024</a>, are removed. Often this leaves only one 'superior' replicate probe in the group.
2. If there is more than one probe left in the group after removing inferior probes*, a mean is taken across the remaining probes. *Reasons for this include that one probe may have the best sensitivity while another has the best precision, group mean may be superior to any one probe, or probes may have no inferior/superior label due to insufficient evidence.

<b>Note this function also removes all probes mapped to "chromosome 0" in the Illumina manifest.</b>
