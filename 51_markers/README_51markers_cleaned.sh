##################################
##################################
#51 clean marker pipeline
##################################
##################################
#Objective: Run phylogenetic analyses using a custom set of reference genomes and using a set of 51 selected marker genes

#General notice: The most important thing is a well named a short but descriptive genome name (DO NOT USE '-' and limit use of other symbols). For the set of reference genomes I chose to  use the GCA_XXX number as point of reference and I add it into the fasta header in step 3.

#To add a more helpful name including the taxonomy I have a mapping file (current reference set for the archaea is v5 and for bacteria is v5)



#Version programs
##################################
#prokka 1.14-dev
#Python 2.7.5
#perl v5.16.3
#HMMER 3.1b2
#blastp (Protein-Protein BLAST 2.7.1+)
#diamond 0.9.22
#IQ-TREE multicore version 1.6.7
#MAFFT v7.407
#BMGE-1.12








#1. set working dir and get essential files
##################################
cd /export/lv1/user/spang_team/Projects/ALR_DPANN/Phylogeny/51markers




#2. cp prokka protein files of reference genomes with genomes of interest
##################################
mkdir faa

#cp arcref v5, 355 genomes
cp ~/../spang_team/Databases/Archaea_Reference_Set/v5/prokka/*/*faa faa

#cp genomes of interest, i.e. UAP2, 12 genomes
cp ~/../spang_team/Projects/ALR_DPANN/faa/renamed/*faa faa/

ll faa//*faa | wc -l
#total = 401 genomes





#3. modify the header of the faa file so that it can be split at a later point, therefore we add the file name and separate it with the old sequence header with a '-'
##################################
'''
Because of how the script separates the bin ID and the protein ID:
MAKE SURE THAT THERE IS NO '-' in your bin name
if a '-' is in your bin-name the file will be split incorrectly and might give a random error
(or otherwise change the script accordingly!!!)
'''

cd faa

mkdir renamed

for i in *faa; do awk '/>/{sub(">","&"FILENAME"|");sub(/\.faa/,x)}1' $i > renamed/$i; done

#rm single faa anc create concatenated file
rm *faa

cat renamed/*faa > All_Genomes.faa

cd ..






#4. create list of used genomes
# the file 'renaming' will be later used to generate an ID that can be used to split the genome ID from the contig ID
##################################
cd faa/renamed

ls *faa > ../../List_of_genomes

cd ../..

#create a file for renaming and splitting files after hmmsearchTable
sed 's/.faa//g' List_of_genomes > temp

awk -F'\t' -v OFS='\t' '{print $1, $1}' temp > temp1

#add in a unique ID that can used later for splitting the name
awk  -F'\t' -v OFS='\t'  '{gsub("$","|CUT",$2)}1' temp1 > renaming

rm temp*

rm faa/renamed/*faa







#5. search genes and clean up file names and extract proteins sequences
##################################
#cp list of genes for looping
cp /export/lv1/user/spang_team/Databases/52_tested_markers/52_markers_list .

mkdir Hmmersearch

~/../spang_team/Scripts/Hmmer/hmmsearchTable faa/All_Genomes.faa /export/lv1/user/spang_team/Databases/52_tested_markers/52tested_inclTIGR.hmm 40 -E 1e-5 > Hmmersearch/HMMscan_Output_e5

wc -l  Hmmersearch/HMMscan_Output_e5
#289283

#duplicate column1 (for cosmetics and easier searching later)
awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $1}' Hmmersearch/HMMscan_Output_e5 > Hmmersearch/temp


#rename bin names in column 1 to be able to later split headers
while read from to; do
sed -i "s/$from/$to/" Hmmersearch/temp;
done < renaming

#separate the 51 marker genes into indiv. files
mkdir Hmmersearch/OGs

for sample in `cat  52_markers_list`; do grep "$sample" Hmmersearch/temp > Hmmersearch/OGs/${sample}.txt; done

ll Hmmersearch/OGs/*txt | wc -l
#51

#cut after first column to get only binID
mkdir Hmmersearch/split

#clean header
for sample in `cat  52_markers_list`; do awk -F'\t' -v OFS='\t' '{split($1,a,"\|CUT"); print a[1], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19 }' Hmmersearch/OGs/$sample* > Hmmersearch/split/$sample |LC_ALL=C  sort ; done


#search and remove duplicates and keep only the best hit
#notice, the hmmsearchTable script does check for duplicates, however, here we specifically search for gene duplications hmmsearchTable would not pick up because it compares the gene name and here we compare the genome name
mkdir Hmmersearch/deduplicated

for sample in `cat 52_markers_list`; do perl ~/../spang_team/Scripts/Others/best_hmmer.pl Hmmersearch/split/$sample* Hmmersearch/deduplicated/$sample; done > duplicated.txt


#get list of proteins to extract from faa file
mkdir Hmmersearch/protein_list

for sample in `cat  52_markers_list`; do awk -F'\t' -v OFS='\t' '{print $19 }' Hmmersearch/deduplicated/$sample* > Hmmersearch/protein_list/$sample |LC_ALL=C  sort ; done


#extract faa sequences
mkdir Marker_Genes

for sample in `cat  52_markers_list`; do perl ~/../spang_team/Scripts/Others/screen_list_new.pl Hmmersearch/protein_list/$sample faa/All_Genomes.faa keep > Marker_Genes/${sample}.faa; done

#shorten header to be able to concatenate later
mkdir Marker_Genes/renamed

for sample in `cat 52_markers_list`
do
cut -f1 -d "|" Marker_Genes/$sample*>> Marker_Genes/renamed/${sample}.faa
done








#6. Prepare alignment
#Notice, adjust the parallel option (-j) according to the free cpus on ada, take care not to use more than 30% of the available cpus
##################################
#6a. align with mafft_linsi
mkdir Alignment
mkdir Alignment/mafft_linsi

parallel -j8 'i={}; nice -n 10 /export/lv1/user/spang_team/Scripts/Mafft/bin/mafft-linsi --reorder --thread 4 Marker_Genes/renamed/$i* > Alignment/mafft_linsi/${i}.aln ' ::: `cat 52_markers_list`


#6b. Trim using BMGE
mkdir Alignment/BMGE
mkdir Alignment/BMGE/h0.55

parallel -j2 'i={}; nice -n 10 java -jar /opt/biolinux/BMGE-1.12/BMGE.jar -i Alignment/mafft_linsi/$i* -t AA -m BLOSUM30 -h 0.55 -of Alignment/BMGE/h0.55/$i ' :::  `cat 52_markers_list`



#6c. concatenate
mkdir Alignment/concatenated

/export/lv1/user/spang_team/Scripts/catfasta2phyml/catfasta2phyml.pl -f -c Alignment/BMGE/h0.55/* > Alignment/concatenated/52markers_ALR_DPANN_ArcRefv5_v1.faa






#7. run iqtree
##################################
mkdir Phylogeny
mkdir Phylogeny/IQtree
mkdir Phylogeny/IQtree/v1

cp Alignment/concatenated/52markers_ALR_DPANN_ArcRefv5_v1.faa Phylogeny/IQtree/v1

cd Phylogeny/IQtree/v1

 nice -n 10 iqtree -s 52markers_ALR_DPANN_ArcRefv5_v1.faa -m LG+C60+F+R -nt AUTO -bb 1000 -alrt 1000


#rename treefile for easier handling (mapping file available for v5 set)
perl ~/../spang_team/Scripts/Others/Replace_tree_names.pl ~/../spang_team/Databases/Archaea_Reference_Set/names_to_replace_mappingJan2018 52markers_ALR_DPANN_ArcRefv5_v1.faa.treefile > 52markers_ALR_DPANN_ArcRefv5_v1.faa.treefile_renamed
