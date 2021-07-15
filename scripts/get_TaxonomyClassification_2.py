import sys,os
from Bio import SeqIO

GoodIds_f=sys.argv[1]

TaxR='kingdom,phylum,class,order,family,genus,species'
TaxRanks=TaxR.split(',')

NCBI_taxonomy_Path='~/DBs/Taxonomy/NCBI/new_taxdump/'


NCBI_nucl_Acc2TaxId_file='~/DBs/Taxonomy/NCBI/nucl_gb.accession2taxid'

NCBI_nucl_MergedTaxIds_file='~/DBs/Taxonomy/NCBI/merged.dmp'


def get_GoodIds(GoodIds_f):
    L_Ids=[]
    f=open(GoodIds_f,'r')
    for line in f:
        Acc=line.strip().split('\t')[0]
        L_Ids.append(Acc)
    f.close()
    return L_Ids

Accessions=get_GoodIds(GoodIds_f)


def get_acc_taxId(Acc2Taxid_file):
    '''This function gives the taxid for every accession.version'''
    Acc_Taxid_dict={}
    f=open(Acc2Taxid_file,'r')
    f.readline() #skip the header
    for line in f:
        l=line.strip().split()
        AccVersion=l[1]
        taxId=l[2]
        Acc_Taxid_dict[AccVersion]=taxId
    f.close()
    return Acc_Taxid_dict


AccTaxid_temp=get_acc_taxId(NCBI_nucl_Acc2TaxId_file)

def get_merged_Taxids(NCBI_nucl_MergedTaxIds_file):
    Dtxid_old_new={}
    f=open(NCBI_nucl_MergedTaxIds_file, 'r')
    for line in f:
        l=line.split('|')
        old_txid=l[0].strip()
        new_txid=l[1].strip()
        Dtxid_old_new[old_txid]=new_txid
    f.close()
    return Dtxid_old_new

Dtxid_old_new=get_merged_Taxids(NCBI_nucl_MergedTaxIds_file)


def update_AccTxid(AccTaxid_temp, Dtxid_old_new):
    AccTaxid={}
    for acc in AccTaxid_temp:
        if AccTaxid_temp[acc] in Dtxid_old_new:
            AccTaxid[acc]=Dtxid_old_new[AccTaxid_temp[acc]]
        else:
            AccTaxid[acc]=AccTaxid_temp[acc]
    AccTaxid_NoVersion={}
    for ac in AccTaxid:
        ac1=ac.split('.')[0]
        AccTaxid_NoVersion[ac1]=AccTaxid[ac]
    return AccTaxid_NoVersion

AccTaxid=update_AccTxid(AccTaxid_temp, Dtxid_old_new)

def select_Accs_taxIds(Accessions, AccTaxid):
    D_selectedAcc_Taxid={}
    L_AccsNotFound=[]
    for acc in Accessions:
        acc1=acc.split('.')[0]
        try:
            Acc_taxid=AccTaxid[acc1]
            D_selectedAcc_Taxid[acc]=Acc_taxid
        except KeyError:
            L_AccsNotFound.append(acc)
    if len(L_AccsNotFound)>0:
        f_out=open(GoodIds_f.split('.')[0]+'_AccsNotFoundInTaxonomy.txt', 'w')
        f_out.write('\n'.join(L_AccsNotFound)+'\n')
        f_out.close()
    return D_selectedAcc_Taxid, L_AccsNotFound


S=select_Accs_taxIds(Accessions, AccTaxid)

SelectedAccs_taxids=S[0]
AccsNotFound=S[1]



def get_Taxid_lineage():
    TaxidLineageFile=NCBI_taxonomy_Path+'taxidlineage.dmp'
    D_Taxid_lineage={}
    f_taxidlin=open(TaxidLineageFile,'r')
    for line in f_taxidlin:
        l=line.split('\t|\n')[0].split('\t|\t')
        TaxId=l[0].strip()
        #Org_taxid=l[1]
        lineage_Txid=l[1].strip()
        D_Taxid_lineage[TaxId]=lineage_Txid+' '+TaxId
    f_taxidlin.close()
    return D_Taxid_lineage

Taxid_Lineage=get_Taxid_lineage()


def get_TaxId_fullNamesLineage():
    D_Taxid_lineageNames={}
    Taxid_FullNamesLineage=NCBI_taxonomy_Path+'fullnamelineage.dmp'
    f_fl=open(Taxid_FullNamesLineage, 'r')
    for line in f_fl:
        l=line.split('\t|\n')[0].split('\t|\t')
        Taxid_Organism=l[0].strip()
        sciName_org=l[1].strip()
        TaxIds_lineage=l[2].strip()
        D_Taxid_lineageNames[Taxid_Organism]=TaxIds_lineage+sciName_org
    f_fl.close()
    return D_Taxid_lineageNames

Taxid_lineageNames=get_TaxId_fullNamesLineage()


def get_TaxId_ranks():
    D_Taxid_Ranks={}
    taxId_ranksFile=NCBI_taxonomy_Path+'nodes.dmp'
    f_ranks=open(taxId_ranksFile,'r')
    for line in f_ranks:
        l=line.split('\t|\t')
        taxId_r=l[0].strip()
        rank=l[2].strip()
        D_Taxid_Ranks[taxId_r]=rank
    f_ranks.close()
    return D_Taxid_Ranks


Taxid_Ranks=get_TaxId_ranks()



def get_SelectedAccs_rankedLineage(SelectedAccs_taxids, Taxid_lineageNames, Taxid_Lineage, Taxid_Ranks, TaxRanks):
    D_Acc_rankedTaxids={}
    D_Acc_rankedlineage={}
    L_DelNodes=[]
    for acc in SelectedAccs_taxids:
        L_TaxidsLineage=[]
        L_TaxidsLineageNames=[]
        try:
            Acc_taxid=SelectedAccs_taxids[acc]
            tid_lineName=Taxid_lineageNames[Acc_taxid].split(';')
            tid_lineage=Taxid_Lineage[Acc_taxid].split()
            for index,taxid in enumerate(tid_lineage):
                rank=Taxid_Ranks[taxid]
                if rank in TaxRanks:
                    L_TaxidsLineage.append(rank+'_'+taxid)
                    lineageName=tid_lineName[index].strip()
                    L_TaxidsLineageNames.append(rank+'_'+lineageName)
                    D_Acc_rankedTaxids[acc]=L_TaxidsLineage
                    D_Acc_rankedlineage[acc]=L_TaxidsLineageNames
        except KeyError:
            L_DelNodes.append(acc)
    if len(L_DelNodes)>0:
        f_out_checkDelNodes=open(GoodIds_f.split('.')[0]+'_CheckDelNodes.txt','w')
        f_out_checkDelNodes.write('\n'.join(L_DelNodes)+'\n')
        f_out_checkDelNodes.close()
    D_Acc_rankedlineage_Final={}
    D_Acc_rankedTaxids_Final={}
    for Acc in D_Acc_rankedlineage:
        Acc_Lineage=D_Acc_rankedlineage[Acc]
        King=Acc_Lineage[0]
        if King=='kingdom_Metazoa':
            D_Acc_rankedlineage_Final[Acc]=Acc_Lineage
            D_Acc_rankedTaxids_Final[Acc]=D_Acc_rankedTaxids[Acc]        
    return D_Acc_rankedlineage_Final, D_Acc_rankedTaxids_Final



SelectedAccs_rankedLineage=get_SelectedAccs_rankedLineage(SelectedAccs_taxids, Taxid_lineageNames, Taxid_Lineage, Taxid_Ranks, TaxRanks)
E=SelectedAccs_rankedLineage[0]
F=SelectedAccs_rankedLineage[1]

print E[E.keys()[0]]
print F[F.keys()[0]]


def writeTaxonomyFile(GoodIds_f, E, F, TaxRanks):
    f_TaxL=open(GoodIds_f.split('.')[0]+'_Taxonomy_Lineages.txt', 'w')
    f_TaxT=open(GoodIds_f.split('.')[0]+'_Taxonomy_TaxIDs.txt','w')
    for Id in E:
        TaxLineage=E[Id]
        TaxTaxid=F[Id]
        for index, rank in enumerate(TaxRanks):
            try:
                if TaxLineage[index].split('_')[0] != rank:
                    TaxLineage.insert(index, rank+'_NA')
            except IndexError:
                TaxLineage.append(rank+'_NA')
            try:
                if TaxTaxid[index].split('_')[0] != rank:
                    TaxTaxid.insert(index, rank+'_NA')
            except IndexError:
                TaxTaxid.append(rank+'_NA')
        f_TaxL.write(Id+'\t'+ '\t'.join(TaxLineage)+ '\n')
        f_TaxT.write(Id+'\t'+ '\t'.join(TaxTaxid) +'\n')
    f_TaxL.close()
    f_TaxT.close()

writeTaxonomyFile(GoodIds_f, E, F,TaxRanks)






