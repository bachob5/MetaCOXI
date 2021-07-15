#/usr/local/bin python2.7
# requires Linux server multithread (preferably 16) and 16 GB of RAM


####Extract ENA data####

wget http://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/std/*

python writeNuc_fromENAFiles.py

#Translate into amino acids. Requires newTranslator_CExtract.py script
python wrapper_ExecuteNewTraslator.py


chmod +x runHMMsearch.sh
./runHMMsearch.sh

python HMMparser_COXI_Extractor.py

python WritePotentialEntries_NCBI.py

python compare_GB_PFAM_annotations_2.py

python get_COXISynonyms_2.py

python WriteFinalNuc_rmFalsePSynonyms.py
cat *_rmFPsynFinalIDs.txt > AllENAr142_SelectedEntries_FinalIDs.temp.txt

#Detect pseudogenes
python Check_Pseudogenes_NCBIFlatF_2.py

# If you do not find pseudogenes you can add an empty file called Pseugened_Accs.txt to get the subsequent code works
#touch Pseudogenes_Accs.txt

cat *Pseudogenes_Entries.txt > Pseudogenes_Accs.txt

python Remove_PseudoEntries_2.py AllENAr142_SelectedEntries_FinalIDs.temp.txt Pseudogenes_Accs.txt 

python get_TaxonomyClassification_2.py AllENAr142_SelectedEntries_FinalIDs.txt 

cat *_Good_Excluded_Entries_SelectedGoodCOX1_rmFPsynonyms.fasta > AllENAr142_SelectedEntries_FinalIDs.temp.fasta

python writeFinalNucTaxo_withoutNs_afterTaxCheck.py AllENAr142_SelectedEntries_FinalIDs.temp.fasta 100 5 AllENAr142_SelectedEntries_FinalIDs_CheckDelNodes.txt AllENAr142_SelectedEntries_FinalIDs_AccsNotFoundInTaxonomy.txt Pseudogenes_Accs.txt




####Extract BOLD Data####

python writeBoldExtractor.py

#Extract sequences
python writeBoldNucTaxonomy_ByMarker1.py 100 5

#Join BOLD records
cat *COI-3P_COI-5P_joinedRecs.fasta Arthropods/*COI-3P_COI-5P_joinedRecs.fasta > AllBoldRecs_COXI.fasta

python generate_Separate_TaxonomyFasta_Bold.py AllBoldRecs_COXI.fasta

python map_BoldTaxonomyOnNCBITaxonomy_2.py AllBoldRecs_COXI_Taxonomy.tsv

python select_BoldAccsNotinNCBI.py AllBoldRecs_COXI_Taxonomy_AccsNotFoundInTaxonomy.txt AllBoldRecs_COXI_Taxonomy.tsv

python map_Taxonomy_BOLDrecsNotinNCBI_2.py AllBoldRecs_COXI_Taxonomy_AccsNotFoundInTaxonomy_Accs2AnalyzeTaxonomy.txt

python Split_LargeFasta_runTranslator_2.py

./runHMMsearch.sh 

python HMMparser_COXI_Extractor.py

cat *_codeframes_MatchIDs_PfamProfile.txt temp_AllBoldRecs_COXI_Taxonomy_*/*_codeframes_MatchIDs_PfamProfile.txt > BOLD_PassedCOXIPFAM.txt

cat AllBoldRecs_COXI_Taxonomy_Taxonomy_Lineages.txt AllBoldRecs_COXI_Taxonomy_AccsNotFoundInTaxonomy_Accs2AnalyzeTaxonomy_Taxonomy_Lineages.txt > Final_AllBoldRecs_COXI_mappedNCBITaxonomy.txt

python Generate_FinalBOLDFiles_2.py AllBoldRecs_COXI_SHortID.fasta Final_AllBoldRecs_COXI_mappedNCBITaxonomy.txt AllBoldRecs_COXI_Taxonomy_AccsNotFoundInTaxonomy_Accs2AnalyzeTaxonomy.txt FastaOUT/BOLD_PassedCOXIPFAM.txt



####JOIN BOLD & ENA####

python Join_BOLD_ENA_bis_2.py FinalBOLDRecs.fasta AllENAr142_SelectedEntries_FinalIDs_WithoutNs.fasta FinalBOLDRecs_taxonomy.tsv AllENAr142_SelectedEntries_FinalIDs_Taxonomy_Lineages_WithoutNNs.tsv
python Join_TaxIDs_BOLD_ENA.py BOLD_ENA_Dereplicated_COXI_Taxonomy.tsv FinalBOLDRecs_taxonomy_TaxIDs.tsv AllENAr142_SelectedEntries_FinalIDs_Taxonomy_TaxIDs_WithoutNNs.tsv


#CreateUniqueBoldTaxonomy.py
python CreateUniqueBoldTaxonomy_2.py

#Generate Metadata
cat *_Good_Excluded_Entries.gb > All_Good_Excluded_Entries_ENArel142.gb
python getMetadata_ENA.py BOLD_ENA_Dereplicated_COXI_Taxonomy.tsv All_Good_Excluded_Entries_ENArel142.gb

#Generate Metadata BOLD
python getMetadata_BOLD_2.py IDs2SearchBold.tsv All_FullRecords_BOLD.tsv

#Final command
python Join_Metadata_ENA_BOLD.py




