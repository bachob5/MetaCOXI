import sys,os
from Bio import SeqIO

GoodIds_f=sys.argv[1]

TaxR='kingdom,phylum,class,order,family,genus,species'
TaxRanks=TaxR.split(',')

selected_Divs=['1','2','5','6','10']

NCBI_taxonomy_Path='~/DBs/Taxonomy/NCBI/new_taxdump/'

NCBI_nucl_Acc2TaxId_file='~/DBs/Taxonomy/NCBI/nucl_gb.accession2taxid'
NCBI_nucl_MergedTaxIds_file='~/DBs/Taxonomy/NCBI/merged.dmp'

#@jit(nopython=True)
def get_GoodIds(GoodIds_f):
    L_Ids=[]
    f=open(GoodIds_f,'r')
    for line in f:
        Acc=line.strip().split('\t')[0]
        L_Ids.append(Acc)
    f.close()
    return L_Ids

Accessions=get_GoodIds(GoodIds_f)

#@jit(nopython=True)
def get_acc_taxId(Acc2Taxid_file):
    '''This function gives the taxid for every accession.version'''
    Acc_Taxid_dict={}
    f=open(Acc2Taxid_file,'r')
    f.readline() #skip the header
    for line in f:
        l=line.strip().split()
        AccVersion=l[1].split('.')[0]
        taxId=l[2]
        Acc_Taxid_dict[AccVersion]=taxId
    f.close()
    return Acc_Taxid_dict


AccTaxid_temp=get_acc_taxId(NCBI_nucl_Acc2TaxId_file)

#@jit(nopython=True)
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

#@jit(nopython=True)
def update_AccTxid(AccTaxid_temp, Dtxid_old_new):
    AccTaxid={}
    for acc in AccTaxid_temp:
        if AccTaxid_temp[acc] in Dtxid_old_new:
            AccTaxid[acc]=Dtxid_old_new[AccTaxid_temp[acc]]
        else:
            AccTaxid[acc]=AccTaxid_temp[acc]
    return AccTaxid

AccTaxid=update_AccTxid(AccTaxid_temp, Dtxid_old_new)

#@jit(nopython=True)
def select_Accs_taxIds(Accessions, AccTaxid):
    D_selectedAcc_Taxid={}
    L_AccsNotFound=[]
    for acc in Accessions:
        try:
            Acc_taxid=AccTaxid[acc]
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
#print SelectedAccs_taxids

#@jit(nopython=True)
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
print Taxid_Lineage.items()[:20]

#@jit(nopython=True)
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
print Taxid_lineageNames.items()[:20]

#@jit(nopython=True)
def get_TaxId_ranks():
    D_Taxid_Ranks={}
    D_Taxid_MitoGenCode={}
    taxId_ranksFile=NCBI_taxonomy_Path+'nodes.dmp'
    f_ranks=open(taxId_ranksFile,'r')
    for line in f_ranks:
        l=line.split('\t|\t')
        taxId_r=l[0].strip()
        rank=l[2].strip()
        Div=l[4].strip()
        MitoCode=l[8].strip()
        if Div in selected_Divs:
            D_Taxid_Ranks[taxId_r]=rank
            D_Taxid_MitoGenCode[taxId_r]=MitoCode
    f_ranks.close()
    return D_Taxid_Ranks, D_Taxid_MitoGenCode


tx_Rank_MitoCode=get_TaxId_ranks()
Taxid_Ranks=tx_Rank_MitoCode[0]
tx_MitoCode=tx_Rank_MitoCode[1]
#print tx_MitoCode.items()[:90]


#@jit(nopython=True)
def get_SelectedAccs_rankedLineage(SelectedAccs_taxids, Taxid_lineageNames, Taxid_Lineage, Taxid_Ranks, TaxRanks, tx_MitoCode):
    D_Acc_rankedTaxids={}
    D_Acc_rankedlineage={}
    D_Acc_MitoCode={}
    L_DelNodes=[]
    for acc in SelectedAccs_taxids:
        L_TaxidsLineage=[]
        L_TaxidsLineageNames=[]
        try:
            Acc_taxid=SelectedAccs_taxids[acc]
            tid_lineName=Taxid_lineageNames[Acc_taxid].split(';')
            tid_lineage=Taxid_Lineage[Acc_taxid].split()
            for index,taxid in enumerate(tid_lineage):
                try:
                    rank=Taxid_Ranks[taxid]
                    print acc, rank
                    if rank in TaxRanks:
                        L_TaxidsLineage.append(rank+'_'+taxid)
                        lineageName=tid_lineName[index].strip()
                        print lineageName
                        L_TaxidsLineageNames.append(rank+'_'+lineageName)
                        D_Acc_rankedTaxids[acc]=L_TaxidsLineage
                        D_Acc_rankedlineage[acc]=L_TaxidsLineageNames
                except KeyError:
                    print acc
        except KeyError:
            L_DelNodes.append(acc)
    for acc_tx in D_Acc_rankedTaxids:
        lastTaxon_taxid=D_Acc_rankedTaxids[acc_tx][-1].split('_')[-1]
        try:
            mCode=tx_MitoCode[lastTaxon_taxid]
            D_Acc_MitoCode.setdefault(mCode, []).append(acc_tx)
        except KeyError:
            print '1->', acc_tx
        #D_Acc_MitoCode[acc]=tx_MitoCode[L_TaxidsLineage[-1]]
    if len(L_DelNodes)>0:
        f_out_checkDelNodes=open(GoodIds_f.split('.')[0]+'_CheckDelNodes.txt','w')
        f_out_checkDelNodes.write('\n'.join(L_DelNodes)+'\n')
        f_out_checkDelNodes.close()          
    return D_Acc_rankedlineage, D_Acc_rankedTaxids, D_Acc_MitoCode



SelectedAccs_rankedLineage=get_SelectedAccs_rankedLineage(SelectedAccs_taxids, Taxid_lineageNames, Taxid_Lineage, Taxid_Ranks, TaxRanks, tx_MitoCode)
E=SelectedAccs_rankedLineage[0]
F=SelectedAccs_rankedLineage[1]
Acc_MitoCode=SelectedAccs_rankedLineage[2]

#print E
#print F

#@jit(nopython=True)
def writeAcc_Mitocode(MitoCode_Acc):
    for c in MitoCode_Acc:
        if len(MitoCode_Acc[c])>0:
            f_out=open(GoodIds_f.split('.')[0]+'_'+c+'_withMitocode.txt','w')
            f_out.write('\n'.join(list(set(MitoCode_Acc[c])))+'\n')
            f_out.close()
            
Mcode=writeAcc_Mitocode(Acc_MitoCode)


#@jit(nopython=True)
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

