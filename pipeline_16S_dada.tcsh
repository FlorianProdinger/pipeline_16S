#This is the 1st version of a 16S analysis pipeline
#I want to basically copy my 18S pipeline, but for 16S


##########################################################################
#QSUB COMMANDS 
#PBS -q cdb 
#PBS -N 16S_dada_Flo
#PBS -l mem=50gb  
#PBS -l ncpus=8
#PBS -o qsub_out_messages.txt
#PBS -e qsub_error_messages.txt
#PBS -m ea
#PBS -M florian.prodinger@gmx.net
##########################################################################

#Loading necessary modules fot his sript
#qiime2/2020.2
#eval `/usr/bin/modulecmd tcsh load qiime2/2018.11`
eval `/usr/bin/modulecmd tcsh load qiime2/2020.2`
eval `/usr/bin/modulecmd tcsh load Python/3.6.5`

### DEFINING VARIABLES

set OUT_DIR="/lustre1/aptmp/florian/uranouchi/16S_uranouchi/20200630_dada_ASV_90_vsearch/"
#"/lustre1/aptmp/florian/uranouchi/16S_uranouchi/20200525_dada_all_data/"
#"/lustre1/aptmp/florian/uranouchi/16S_uranouchi/20190805_version1_more_data/"
#"/lustre1/aptmp/florian/uranouchi/16S_uranouchi/20190718_version1/"
set MANIFEST="qiime_manifest.txt"
set METADATA_FILE="qiime_metadata.tsv"
#"qiime_metadata.csv"

set CLUST_PER = "97"
set NPCUS_THREAD="8"


set SILVA_SEQ_DIR = "/lustre1/aptmp/florian/uranouchi/16S_uranouchi/SILVA/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/{$CLUST_PER}/"
set SIVLA_TAX_DIR = "/lustre1/aptmp/florian/uranouchi/16S_uranouchi/SILVA/SILVA_132_QIIME_release/taxonomy/taxonomy_all/{$CLUST_PER}/"
set SILVA_FNA = "silva_132_{$CLUST_PER}_16S.fna"
set SILVA_TAX = "majority_taxonomy_all_levels.txt"

##The MANIFEST file can be gerenated with a script, but manuall creation is recommended.
## tutorial: https://docs.qiime2.org/2018.11/tutorials/importing/
##The METADATA_FILE file is needed for a barplot and other plots and has to be created manually
## tutorial: https://docs.qiime2.org/2018.11/tutorials/metadata/
###########################################################################


##  1  IMPORT
This script imports the files specified in the "manifest" 
echo "import using $MANIFEST"
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ${OUT_DIR}$MANIFEST \
--output-path ${OUT_DIR}out1-import.qza \
--input-format PairedEndFastqManifestPhred33


#### DADA2 includes:
#  PRIMER REMOVAL
#  MERGING
#  QUALITY FILTER
#  DEREPLICATION
#  CHIMERA CHECK

#primers used
#set PRIMER_F="CCTACGGGNBGCASCAG"
#set PRIMER_R="GACTACNVGGGTATCTAATCC"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ${OUT_DIR}out1-import.qza \
  --p-trunc-len-f 250\
  --p-trunc-len-r 250\
  --p-trim-left-f 17  \
  --p-trim-left-r 22  \
  --p-n-threads $NPCUS_THREAD \
  --o-denoising-stats ${OUT_DIR}out_dada2_denoising-stats \
  --o-table ${OUT_DIR}out5-derep-table.qza \
  --o-representative-sequences ${OUT_DIR}out5-derep-sequences.qza

################# SILVA IMPORT ################
###
echo "SILVA sequence import at $CLUST_PER%"
##6.1) importing 97% sequence data
qiime tools import --type 'FeatureData[Sequence]' \
--input-path $SILVA_SEQ_DIR$SILVA_FNA  \
--output-path ${OUT_DIR}out6.1-silva_sequence_0.{$CLUST_PER}_16Sonly.qza
##
echo "SILVA $SILVA_TAX import ($CLUST_PER%)"
##6.2) import 97% taxonomy data
qiime tools import  --type 'FeatureData[Taxonomy]' \
--input-path $SIVLA_TAX_DIR$SILVA_TAX \
--input-format HeaderlessTSVTaxonomyFormat \
--output-path ${OUT_DIR}out6.2-silva_0.{$CLUST_PER}_taxonomy_majority_out.qza
################# SILVA IMPORTED ################



# 8.3) tree generation
#   ALIGNMENT & TREE
#--i-sequences ${OUT_DIR}out5-derep-table.qza \
echo "alignment & tree generation"
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences ${OUT_DIR}out5-derep-sequences.qza \
--o-alignment ${OUT_DIR}out8-seqs-nonchimeras_clustered_{$CLUST_PER}-aligned.qza \
--o-masked-alignment ${OUT_DIR}out8-seqs-nonchimeras_clustered_{$CLUST_PER}-masked-aligned.qza \
--o-tree ${OUT_DIR}out8-seqs-nonchimeras_clustered_{$CLUST_PER}_unrooted-tree.qza \
--o-rooted-tree ${OUT_DIR}out8-seqs-nonchimeras_clustered_{$CLUST_PER}_rooted-tree.qza

##  9  TAXONOMIC ANOTATION
#This commands annotates reads
#--i-query ${OUT_DIR}out5-derep-table.qza \
echo "vsearch feature classifier at $CLUST_PER%"
qiime feature-classifier classify-consensus-vsearch \
--i-query ${OUT_DIR}out5-derep-sequences.qza \
--i-reference-taxonomy ${OUT_DIR}out6.2-silva_0.{$CLUST_PER}_taxonomy_majority_out.qza \
--i-reference-reads ${OUT_DIR}out6.1-silva_sequence_0.{$CLUST_PER}_16Sonly.qza \
--o-classification ${OUT_DIR}out9-classify-vsearch_{$CLUST_PER}p.qza \
--p-perc-identity 0.9 \
--p-maxaccepts 1 \
--p-threads $NPCUS_THREAD

 
#  9.1 TAXONOMIC OTU FILTER
echo "filtering singeltons"
#${OUT_DIR}out8-derep-filtered_table_clustered_{$CLUST_PER}.qza\
qiime feature-table filter-features \
--i-table ${OUT_DIR}out5-derep-table.qza\
--m-metadata-file ${OUT_DIR}out9-classify-vsearch_{$CLUST_PER}p.qza \
--p-where "Taxon NOT LIKE '%Unassigned%'"\
--p-min-frequency 2 \
--o-filtered-table ${OUT_DIR}out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton.qza


##    a) ASV - table  
qiime feature-table summarize \
--i-table ${OUT_DIR}out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton.qza\
--o-visualization ${OUT_DIR}out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton_OTU.qzv\
--m-sample-metadata-file ${OUT_DIR}$METADATA_FILE
##    b) ASV - sequence list
qiime feature-table tabulate-seqs \
--i-data ${OUT_DIR}out5-derep-sequences.qza \
--o-visualization ${OUT_DIR}out8.4-out5-derep-sequences_{$CLUST_PER}_ASV_Seqs.qzv 


#    e) taxonomic read table as tsv
mkdir -p ${OUT_DIR}export_for_OTU_table
foreach LEVEL (2 3 4 5 6 7 8 9 10)
 qiime taxa collapse --i-table ${OUT_DIR}out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton.qza \
 --i-taxonomy ${OUT_DIR}out9-classify-vsearch_{$CLUST_PER}p.qza\
 --p-level $LEVEL\
 --o-collapsed-table ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton_taxonomy_lvl{$LEVEL}.qza
 qiime tools export --input-path ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton_taxonomy_lvl{$LEVEL}.qza\
 --output-path ${OUT_DIR}export_for_OTU_table
 /bin/mv ${OUT_DIR}export_for_OTU_table/feature-table.biom ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton_taxonomy_lvl{$LEVEL}_feature-table.biom
 biom convert -i ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton_taxonomy_lvl{$LEVEL}_feature-table.biom\
 -o ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton_taxonomy_lvl{$LEVEL}.tsv\
 --to-tsv
end

#    f) OTU table with reads as tsv
mkdir -p ${OUT_DIR}export_for_OTU_table
qiime tools export --input-path ${OUT_DIR}out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton.qza \
--output-path ${OUT_DIR}export_for_OTU_table
/bin/mv ${OUT_DIR}export_for_OTU_table/feature-table.biom ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton_feature-table.biom
biom convert -i ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton_feature-table.biom \
-o ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton_OTUs.tsv \
--to-tsv




##    g) OTU table with reads as tsv and taxonomic annotation


echo "exporting classification (OTU - tax list) from qiime..."
qiime tools export \
 --input-path ${OUT_DIR}out9-classify-vsearch_{$CLUST_PER}p.qza \
 --output-path ${OUT_DIR}export_tax_OTU_list
echo "exporting OTU table without annotation from qiime..."
qiime tools export \
 --input-path  ${OUT_DIR}out5-derep-table.qza \
 --output-path ${OUT_DIR}export_tax_OTU_list

mkdir -p ${OUT_DIR}export_tax_OTU_list

echo "using 'command_rename_first_line_taxonomy_file.R' to rename the classification (OTU - tax) list header"
#replace fist line of taxonomy.tsv with "#OTUID       taxonomy        confidence"
Rscript /lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/pipeline_commands/command_rename_first_line_taxonomy_file.R \
${OUT_DIR}export_tax_OTU_list/taxonomy.tsv \
{$CLUST_PER}

echo "adding metadata (taxonomy) to the OTU table (feature-table.biom)"
biom add-metadata \
 -i ${OUT_DIR}export_tax_OTU_list/feature-table.biom \
 --observation-metadata-fp ${OUT_DIR}export_tax_OTU_list/taxonomy.tsv${CLUST_PER}_renamed.tsv \
 --sc-separated taxonomy \
 -o ${OUT_DIR}export_tax_OTU_list/table-with-taxonomy${CLUST_PER}.biom

echo "exporting the OTU table (feature-table.biom) WITH the added caddification (tax)"
biom convert \
 -i ${OUT_DIR}export_tax_OTU_list/table-with-taxonomy${CLUST_PER}.biom \
 --to-tsv \
 --header-key taxonomy \
 -o ${OUT_DIR}export_tax_OTU_list/out5-derep-table_${CLUST_PER}_with_tax.tsv
/bin/rm -r -f ${OUT_DIR}export_tax_OTU_list/taxonomy.tsv
/bin/rm -r -f ${OUT_DIR}export_tax_OTU_list/feature-table.biom

echo "file: ${OUT_DIR}export_tax_OTU_list/out8-derep-filtered_table_clustered_${CLUST_PER}_with_tax.tsv"


#creating the barplot
echo "creating barplot - PLOTTING"
qiime taxa barplot \
--i-table ${OUT_DIR}out9-derep-clustered_${CLUST_PER}-final_table_no_singelton.qza \
--i-taxonomy ${OUT_DIR}out9-classify-vsearch_${CLUST_PER}p.qza \
--m-metadata-file ${OUT_DIR}$METADATA_FILE \
--o-visualization ${OUT_DIR}out9_classify-table_barplot_ASV.qzv


foreach SUBSAMPLING (1000 2500 5000 50000)
   echo "creating rare faction with $SUBSAMPLING subsampled reads of the dataset clustered at $CLUST_PER%..."
   qiime diversity alpha-rarefaction \
   --i-table ${OUT_DIR}out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton.qza\
   --p-max-depth {$SUBSAMPLING} \
   --m-metadata-file ${OUT_DIR}$METADATA_FILE\
   --o-visualization ${OUT_DIR}out9-derep-culstered_{$CLUST_PER}-final_${SUBSAMPLING}.qzv \
   --p-steps 60

   echo "creating MDS type data with $SUBSAMPLING subsampled reads of the $CLUST_PER dataset..."
   set OUT_FOLDER_NAME="core-metrics-phylogenetic-output_${SUBSAMPLING}_depth_${CLUST_PER}"
   /bin/rm -r -f $OUT_FOLDER_NAME
   qiime diversity core-metrics-phylogenetic \
   --i-table ${OUT_DIR}out9-derep-clustered_{$CLUST_PER}-final_table_no_singelton.qza \
   --i-phylogeny ${OUT_DIR}out8-seqs-nonchimeras_clustered_{$CLUST_PER}_rooted-tree.qza\
   --p-sampling-depth $SUBSAMPLING \
   --m-metadata-file ${OUT_DIR}$METADATA_FILE \
   --output-dir ${OUT_DIR}$OUT_FOLDER_NAME
   ####--p-n-jobs 16 \ #this does not work > "segmentation fault" error
end

