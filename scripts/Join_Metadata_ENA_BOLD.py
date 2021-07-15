from Bio import SeqIO

PathToMetadata='~/Joined_BOLD_ENAr142/Metadata/'

metadata_ENA= PathToMetadata+'Metadata_ENA.tsv'
metadata_BOLD= PathToMetadata+'Metadata_BOLD.tsv'

GeneticCodes_ENA=PathToMetadata+'ENA_IDs_GeneticCodes.tsv'
GeneticCodes_BOLD=PathToMetadata+'BOLD_IDs_GeneticCodes.tsv'

source_ENA='AllENAr142_SelectedEntries_FinalIDs_Taxonomy_Lineages_WithoutNNs.tsv'
source_BOLD='FinalBOLDRecs_taxonomy.tsv'

Final_coreDB_Seqs='BOLD_ENA_Dereplicated_COXI_Seqs.fasta'
Final_coreDB_Taxonomy='BOLD_ENA_Dereplicated_COXI_Taxonomy.tsv'
Final_coreDB_Taxids='BOLD_ENA_Dereplicated_COXI_Taxonomy_Taxids.tsv'


def get_MitoGC(f_GeneticCodes):
    D_Acc_MitoGC={}
    f=open(f_GeneticCodes, 'r')
    #f.readline()
    for line in f:
        l=line.strip().split('\t')
        D_Acc_MitoGC[l[0]]=l[1]
    f.close()
    return D_Acc_MitoGC

BOLD_MitoGC=get_MitoGC(GeneticCodes_BOLD)
ENA_MitoGC=get_MitoGC(GeneticCodes_ENA)

def getMetadata(f_Metadata):
    D_ID_Metadata={}
    D_Acc_ID={}
    f=open(f_Metadata,'r')
    f.readline()
    for line in f:
        l=line.strip().split('\t')
        Acc=l[0]
        Acc_temp=Acc.split(';')
        aMetadata='\t'.join(l[1:])
        if len(Acc_temp)==1:
            D_ID_Metadata[Acc]=aMetadata
        elif len(Acc_temp)==2:
            D_Acc_ID[Acc_temp[1]]=Acc_temp[0]
            D_ID_Metadata[Acc_temp[0]]=aMetadata
    f.close()
    return D_ID_Metadata, D_Acc_ID

ENA_Metadata=getMetadata(metadata_ENA)[0]
Bold_Metadata_temp=getMetadata(metadata_BOLD)
Bold_Metadata=Bold_Metadata_temp[0]
Acc_ProcessID_BOLD=Bold_Metadata_temp[1]

def get_SourceEntries(f_Source):
    D_Acc_True={}
    f=open(f_Source,'r')
    for line in f:
        l=line.strip().split('\t')
        Acc=l[0]
        Acc1=l[0].split('.')[0]
        D_Acc_True[Acc]=True
        D_Acc_True[Acc1]=True
    f.close()
    return D_Acc_True

D_ENA_source=get_SourceEntries(source_ENA)
D_BOLD_source=get_SourceEntries(source_BOLD)


def get_Acc_Seq(Final_coreDB_Seqs):
    D_Acc_Seq={}
    f_Fas=SeqIO.parse(Final_coreDB_Seqs, 'fasta')
    for rec in f_Fas:
        D_Acc_Seq[rec.id]=str(rec.seq)
    return D_Acc_Seq

All_Acc_Seqs=get_Acc_Seq(Final_coreDB_Seqs)

def get_Acc_Taxonomy(Final_coreDB_Taxonomy):
    D_Acc_Taxonomy={}
    f=open(Final_coreDB_Taxonomy, 'r')
    for line in f:
        l=line.strip().split('\t')
        Acc=l[0]
        Taxonomy='\t'.join(l[1:])
        D_Acc_Taxonomy[Acc]=Taxonomy
    f.close()
    return D_Acc_Taxonomy

All_Acc_Taxonomy=get_Acc_Taxonomy(Final_coreDB_Taxonomy)
All_Acc_Taxids=get_Acc_Taxonomy(Final_coreDB_Taxids)

def Join_Taxonomy_Taxids(All_Acc_Taxonomy, All_Acc_Taxids):
    D_Acc_TaxoTaxids={}
    for ac in All_Acc_Taxonomy:
        print ac
        L=[]
        Taxo=All_Acc_Taxonomy[ac].split('\t')
        Taxids=All_Acc_Taxids[ac].split('\t')
        for i,j in enumerate(Taxo):
            L.append(Taxo[i])
            L.append(Taxids[i])
        D_Acc_TaxoTaxids[ac]='\t'.join(L)
    return D_Acc_TaxoTaxids

Joined_TaxoTaxids=Join_Taxonomy_Taxids(All_Acc_Taxonomy, All_Acc_Taxids)

D_Id_Seqs={}
D_Id_Metadata={}
for Id in All_Acc_Seqs:
    L=[]
    Taxonomy=Joined_TaxoTaxids[Id]
    L.append(Taxonomy)
    if ENA_MitoGC.has_key(Id):
        L.append(ENA_MitoGC[Id])
    else:
        L.append(BOLD_MitoGC[Id])
    if ENA_Metadata.has_key(Id):
        L.insert(0, Id)
        L.append(ENA_Metadata[Id])
        D_Id_Seqs[Id]=All_Acc_Seqs[Id]
    elif Bold_Metadata.has_key(Id):
        L.insert(0, Id)
        L.append(Bold_Metadata[Id])
        D_Id_Seqs[Id]=All_Acc_Seqs[Id]
    else:
        differentID=Acc_ProcessID_BOLD[Id]
        D_Id_Seqs[differentID]=All_Acc_Seqs[Id]
        L.insert(0, differentID)
        if ENA_Metadata.has_key(differentID):
            L.append(ENA_Metadata[differentID])
        else:
            L.append(Bold_Metadata[differentID])
    if D_ENA_source.has_key(Id) and D_BOLD_source.has_key(Id):
        L.append('ENA&BOLD')
    elif D_ENA_source.has_key(Id) and Id not in D_BOLD_source:
        L.append('ENA')
    elif D_BOLD_source.has_key(Id) and Id not in D_ENA_source:
        L.append('BOLD')
    D_Id_Metadata[L[0]]='\t'.join(L[1:])

f_out_seqs=open('FinalCoreDB_COXI_Seqs_2.fasta', 'w')
f_out_Metadata=open('FinalCoreDB_COXI_Metadata_2.tsv','w')
f_out_Metadata.write('Accession\tSequence_Length\tKingdom\tKingdom_Taxid\tPhylum\tPhylum_Taxid\tClass\tClass_Taxid\tOrder\tOrder_Taxid\tFamily\tFamily_Taxid\tGenus\tGenus_Taxid\tSpecies\tSpecies_Taxid\tMitochondrialGeneticCode\tcountry\tcountry_address\thost\tcollection_year\tcollection_complete_date\tgeographical_coordinates\tRetrievalLink\tProvenance_DB\n')
for ID in D_Id_Seqs:
    f_out_seqs.write('>'+ID+'\n'+D_Id_Seqs[ID]+'\n')
f_out_seqs.close()

for IDm in D_Id_Metadata:
    f_out_Metadata.write(IDm+'\t'+ str(len(D_Id_Seqs[IDm]))+'\t' +D_Id_Metadata[IDm]+'\n')
f_out_Metadata.close()
            
    
        






