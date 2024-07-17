import pandas as pd
import gzip

#names of the dataframes
# just uncomment the pair that need to be processed.
#4 different motifs are available - MA0808.1, MA1121.1, MA0090.3, MA0809.2

#table_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA0808.1_binding_sites.csv"
#table_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA1121.1_binding_sites.csv"
#table_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA0090.3_binding_sites.csv"
#table_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA0809.2_binding_sites.csv"
#table_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA0099.3_binding_sites.csv"
table_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA1128.1_binding_sites.csv"


#result_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA0808.1_binding_sites_results.csv"
#result_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA1121.1_binding_sites_results.csv"
#result_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA0090.3_binding_sites_results.csv"
#result_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA0809.2_binding_sites_results.csv"
#result_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA0099.3_binding_sites_results.csv"
result_address = "/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/results/integration_of_HIC/MA1128.1_binding_sites_results.csv"



df = pd.read_csv(table_address, sep=";")
print(df.head())

# initiate a dataframe to store the overlapping intervals
result_df=pd.DataFrame(columns=["Chr", "Target_Start", "Target_End", "HIC_Reference_Start", "HIC_Reference_End", "HIC_bin2", "HIC_distfoldchange"])
result_df.to_csv(result_address, sep='\t', index=False, mode='w')

for chrname in df.seqnames.unique():
    #chrname = "chr1"
    print(chrname)
    reference_table_path = "".join(["/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/data_public/Hi-C/full_data/MSC/MSC.", chrname, ".distnorm.scale2.gz"])
    with gzip.open(reference_table_path, 'rt') as file:
        reference_df = pd.read_csv(file, sep="\t")
    
    reference_df = reference_df[reference_df["dist_foldchange"] > 2]
    target_df = df[df.seqnames == chrname]

    overlapping_intervals = pd.DataFrame({"Chr":[], "Target_Start":[], "Target_End":[], "HIC_Reference_Start":[], "HIC_Reference_End":[], "HIC_bin2":[], "HIC_distfoldchange":[]})
    for target_start, target_end in zip(target_df.start, target_df.end):
        # print(chrname, target_start, target_end)
        list_to_append = []
        for reference_start, reference_bin2, reference_dist_foldchange in zip(reference_df.bin1, reference_df.bin2, reference_df.dist_foldchange):
            #print(chrname, reference_start, reference_bin2, reference_dist_foldchange)
            if target_end >= reference_start and target_start <= reference_start+5000:
                list_to_append.append([chrname, target_start, target_end, reference_start, reference_start+5000, reference_bin2,  reference_dist_foldchange])
        
        overlapping_intervals = pd.concat([overlapping_intervals, pd.DataFrame(list_to_append, columns=["Chr", "Target_Start", "Target_End", "HIC_Reference_Start", "HIC_Reference_End", "HIC_bin2", "HIC_distfoldchange"])])
    overlapping_intervals.to_csv(result_address, sep='\t', index=False, mode='a', header=False)






