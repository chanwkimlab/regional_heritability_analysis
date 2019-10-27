base_path='data/'
gwas_path='/data01/ch6845/UKB_gwas_neale/'


#ukb_sub_bim_path=base_path+'ukb_imp_plink/ukb_imp_chr{}_v3.bim'
bim_path=base_path+'1000G_plink_EUR/EUR_phase3_chr{}.bim'
bfile_path=base_path+'1000G_plink_EUR/EUR_phase3_chr{}'

"""
phenotypes_both_sexes_file_path=gwas_path+'phenotypes.both_sexes.tsv.gz'
phenotypes_male_file_path=gwas_path+'phenotypes.male.tsv.gz'
phenotypes_female_file_path=gwas_path+'phenotypes.female.tsv.gz'

biomarkers_both_sexes_file_path=gwas_path+'biomarkers.both_sexes.tsv.bgz'
biomarkers_male_file_path=gwas_path+'biomarkers.male.tsv.bgz'
biomarkers_female_file_path=gwas_path+'biomarkers.female.tsv.bgz'
"""
phenotypes_both_sexes_v2_file_path=gwas_path+'phenotypes.both_sexes.v2.tsv.bgz'
phenotypes_male_v2_file_path=gwas_path+'phenotypes.male.v2.tsv.bgz'
phenotypes_female_v2_file_path=gwas_path+'phenotypes.female.v2.tsv.bgz'

variants_file_path=gwas_path+'variants.tsv.bgz'


h2_path=base_path+'ukbb_all_h2univar_results.txt'
h2_par_path=base_path+'ukbb_all_h2part_results.txt'
h2_v2_path=base_path+'ukb31063_h2_topline.02Oct2019.tsv'
variants_wo_cm_path=base_path+'variants_wo_cm.bim'
variants_w_cm_path=base_path+'variants_w_cm.bim'
filter_index_path=base_path+'filter.index'
filter_rsid_path=base_path+'filter.rsid'

ukbb_table_path=base_path+'UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Manifest 201807.csv'
ukbb_table_filtered_path=base_path+'UKBB_sumstats_filtered.csv'

phenotypes_filtered_path=base_path+'phenotypes_filtered.csv'

phenotypes_uni_filtered_path=base_path+'phenotypes_uni_filtered.csv'
phenotypes_par_filtered_path=base_path+'phenotypes_par_filtered.csv'
h2_total_par_filtered_path=base_path+'h2_total_par_filtered.csv'

"""
phenotypes_par_filtered_description_dict_path=base_path+'phenotypes_par_filtered_description_dict.tsv'
pleiotropic_loci_description_dict_path=base_path+'pleiotropic_loci_description_dict.tsv'
correlation_description_dict_path=base_path+'correlation_description_dict.tsv'
"""
description_dict_merge_path=base_path+'description_dict_merge.tsv'

annot_path=base_path+'out_annot/{}.{}.annot'
#premunge_path=base_path+'out_sumstats/{}.premunge'
premunge_path='/home/ch6845/premunge/{}.premunge'
print_snps_path=base_path+'hapmap3_snps/hm.{}.snp'

snplist_path=base_path+'w_hm3.snplist'
sumstats_path=base_path+'out_sumstats/{}'

#LDscore regression
wld_path=base_path+'1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.'
ldsc_path=base_path+'out_final/{}.ldsc'


#estimating LD score
ld_path=base_path+'out_annot/{}.{}'