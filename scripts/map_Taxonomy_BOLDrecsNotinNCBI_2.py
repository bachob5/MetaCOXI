import sys,os
from Bio import SeqIO

GoodIds_f=sys.argv[1]


TaxR='kingdom,phylum,class,order,family,genus,species'
TaxRanks=TaxR.split(',')



def taxonomyLevelsNum():
    D_TaxonomyLevelsNum={}
    TaxR='kingdom,phylum,class,order,family,genus,species'
    TaxRanks=TaxR.split(',')
    for i,n in enumerate(TaxRanks):
        D_TaxonomyLevelsNum[i]=n
    return D_TaxonomyLevelsNum

taxonomyLevelsNum=taxonomyLevelsNum()

print taxonomyLevelsNum

selected_Divs=['1','2','5','6','10']


NCBI_taxonomy_Path='~/DBs/Taxonomy/NCBI/new_taxdump/'


NCBI_nucl_Acc2TaxId_file='~/DBs/Taxonomy/NCBI/nucl_gb.accession2taxid' 

NCBI_nucl_MergedTaxIds_file='~/DBs/Taxonomy/NCBI/merged.dmp'


def get_BOLDTaxNames(GoodIds_f):
    D_Name_Acc={}
    D_Acc_TaxLevelNum={}
    D_Acc_Name4MitoCode={}
    D_Acc_AllTaxonomyNoNA={}
    f=open(GoodIds_f, 'r')
    for line in f:
        l=line.strip().split('\t')
        i=-1
        while l[i]=='NA':
            l.remove(l[i])
        taxN=l[i].replace('_', ' ')
        D_Name_Acc[l[0]]=taxN
        D_Acc_AllTaxonomyNoNA[l[0]]=l[1:]
        if l[1]!='NA':
            D_Acc_Name4MitoCode[l[0]]=l[1]
        else:
            D_Acc_Name4MitoCode[l[0]]=l[2]
        for i,n in enumerate(l[1:]):
            D_Acc_TaxLevelNum[l[0]]=i
    f.close()
    return D_Name_Acc, D_Acc_TaxLevelNum, D_Acc_Name4MitoCode, D_Acc_AllTaxonomyNoNA


Res_BOLDtax=get_BOLDTaxNames(GoodIds_f)
lastTaxon=Res_BOLDtax[0]
lastTaxonIndex=Res_BOLDtax[1]
BOLD_Acc_Name4MitoCode=Res_BOLDtax[2]
BOLD_Acc_AllTaxonomyNoNA=Res_BOLDtax[3]

#print lastTaxon
#print lastTaxonIndex

def get_TaxId_ranks():
    D_Taxid_Ranks={}
    D_Taxid_MitoGenCode={} #8
    D_Taxid_Div={}#4
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
            D_Taxid_Div[taxId_r]=Div
            D_Taxid_MitoGenCode[taxId_r]=MitoCode
    f_ranks.close()
    return D_Taxid_Ranks, D_Taxid_Div, D_Taxid_MitoGenCode


D_all=get_TaxId_ranks()

Taxid_Ranks=D_all[0]
Taxid_Div=D_all[1]
Taxid_MitoCode=D_all[2]
print Taxid_MitoCode.items()[:20]

def getTaxID_Names(NCBITaxonomyPath):
    Dict_ScientificNames_TaxID={}
    Dict_Synonyms_TaxID={}
    f_TaxNames=open(NCBITaxonomyPath+'names.dmp','r')
    for line in f_TaxNames:
        l=line.split('\t')
        tax_id=l[0]
        name_txt=l[2]
        #uniqueName=l[4]
        className=l[6]
        if className =='scientific name':
            Dict_ScientificNames_TaxID.setdefault(name_txt,[]).append(tax_id)
        elif className=='synonym':
            Dict_Synonyms_TaxID.setdefault(name_txt,[]).append(tax_id)
    f_TaxNames.close()
    return Dict_ScientificNames_TaxID, Dict_Synonyms_TaxID 


tx_n_Id=getTaxID_Names(NCBI_taxonomy_Path)
SciName_TaxId=tx_n_Id[0]
Synonyms_TaxId=tx_n_Id[1]

       

def mapNames_NCBI_TaxId(Taxid_Ranks, lastTaxon, lastTaxonIndex, SciName_TaxId,taxonomyLevelsNum, BOLD_Acc_Name4MitoCode):
    D_Id_UniqueTaxID={}
    SciNamesNotFound=[]
    D_Id_HigherTxId={}
    for Id in lastTaxon:
        Taxon=lastTaxon[Id]
        #print Taxon
        #try:
        if SciName_TaxId.has_key(Taxon):
            NcbiTaxId_temp=SciName_TaxId[Taxon]
            if len(NcbiTaxId_temp)==1:
                NcbiTaxId=NcbiTaxId_temp[0]
                TaxLevel=lastTaxonIndex[Id]
                #print TaxLevel
                D_Id_UniqueTaxID.setdefault(Id,[]).append(NcbiTaxId)
            elif len(NcbiTaxId_temp)>1:
                for txid in NcbiTaxId_temp:
                    taxLev=lastTaxonIndex[Id]
                    if Taxid_Ranks.has_key(txid):
                        if Taxid_Ranks[txid]==taxonomyLevelsNum[taxLev]:
                            D_Id_UniqueTaxID.setdefault(Id,[]).append(txid)
        elif Synonyms_TaxId.has_key(Taxon):
            NcbiTaxId_temp=Synonyms_TaxId[Taxon]
            if len(NcbiTaxId_temp)==1:
                NcbiTaxId=NcbiTaxId_temp[0]
                TaxLevel=lastTaxonIndex[Id]
                #print TaxLevel
                D_Id_UniqueTaxID.setdefault(Id,[]).append(NcbiTaxId)
            elif len(NcbiTaxId_temp)>1:
                for txid in NcbiTaxId_temp:
                    taxLev=lastTaxonIndex[Id]
                    if Taxid_Ranks.has_key(txid):
                        if Taxid_Ranks[txid]==taxonomyLevelsNum[taxLev]:
                            D_Id_UniqueTaxID.setdefault(Id,[]).append(txid)
        else:
            SciNamesNotFound.append(Id)
        #except KeyError:
            #SciNamesNotFound.append(Id)
            try:
                HigherNameTxID=SciName_TaxId[BOLD_Acc_Name4MitoCode[Id]]
                D_Id_HigherTxId[Id]=HigherNameTxID
            except KeyError:
                print Id
    if len(SciNamesNotFound)>0:
        f_out=open(GoodIds_f.split('.')[0]+'_ScientificNameNotFound.txt','w')
        f_out.write('\n'.join(SciNamesNotFound))
        f_out.close()
    return D_Id_UniqueTaxID, SciNamesNotFound, D_Id_HigherTxId

Names_TaxIds=mapNames_NCBI_TaxId(Taxid_Ranks, lastTaxon, lastTaxonIndex, SciName_TaxId, taxonomyLevelsNum, BOLD_Acc_Name4MitoCode)
Id_UniqueTaxID=Names_TaxIds[0]
SciNamesNotFound=Names_TaxIds[1]
#print Id_UniqueTaxID
SciNamesNotFound_HigherTaxID=Names_TaxIds[2]

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


def compareTaxonomy_SameTaxName(l, d, l1):
    '''
    l=['1','2'] #taxids having the same name at the same last taxonomic level
    d={'1':[1,2,3,4], '2':[1,2,3,7]} #lineages
    l1=[1,2,3,4] #BoldTaxonomy
    '''
    d1={}
    for i in d:
        n=0
        t=d[i]
        for j in l1:
            if j in t:
                n+=1
        d1[i]=n
    l2=[[y,x] for x,y in d1.items()]
    l2.sort(reverse=True)
    return l2[0][1]

def get_SelectedAccs_rankedLineage(SelectedAccs_taxids, Taxid_lineageNames, Taxid_Lineage, Taxid_Ranks, TaxRanks, Taxid_MitoCode, SciNamesNotFound_HigherTaxID, BOLD_Acc_AllTaxonomyNoNA):
    D_Acc_rankedTaxids={}
    D_Acc_rankedlineage={}
    L_DelNodes=[]
    D_MitoCode_Acc={}
    TaxaWithSameNames=[]
    for acc in SelectedAccs_taxids:
        #print acc, SelectedAccs_taxids[acc]
        L_TaxidsLineage=[]
        L_TaxidsLineageNames=[]
        try:
            Acc_taxid_temp=SelectedAccs_taxids[acc]
            #print 'Acc_taxid', acc, Acc_taxid
            if len(Acc_taxid_temp)==1:
                Acc_taxid=Acc_taxid_temp[0]
                tid_lineName=Taxid_lineageNames[Acc_taxid].split(';')
                #print 'tid_lineName', tid_lineName
                tid_lineage=Taxid_Lineage[Acc_taxid].split()
                #print tid_lineage
                for index,taxid in enumerate(tid_lineage):
                    try:
                        rank=Taxid_Ranks[taxid]
                        #print acc, rank
                        if rank in TaxRanks:
                            L_TaxidsLineage.append(rank+'_'+taxid)
                            #print L_TaxidsLineage
                            lineageName=tid_lineName[index].strip()
                            L_TaxidsLineageNames.append(rank+'_'+lineageName)
                            D_Acc_rankedTaxids[acc]=L_TaxidsLineage
                            D_Acc_rankedlineage[acc]=L_TaxidsLineageNames
                    except KeyError:
                        print 'acc, 1txID', acc
            else:
                BoldTaxonomy=BOLD_Acc_AllTaxonomyNoNA[acc]
                d1={}
                for acc_tx in Acc_taxid_temp:
                    tx_lineage=Taxid_lineageNames[acc_tx].split(';')
                    d1[acc_tx]=tx_lineage
                Good_TaxId=compareTaxonomy_SameTaxName(Acc_taxid_temp, d1, BoldTaxonomy)
                tid_lineName=Taxid_lineageNames[Good_TaxId].split(';')
                tid_lineage=Taxid_Lineage[Good_TaxId].split()
                for index,taxid in enumerate(tid_lineage):
                    try:
                        rank=Taxid_Ranks[taxid]
                        #print acc, rank
                        if rank in TaxRanks:
                            L_TaxidsLineage.append(rank+'_'+taxid)
                            #print L_TaxidsLineage
                            lineageName=tid_lineName[index].strip()
                            L_TaxidsLineageNames.append(rank+'_'+lineageName)
                            D_Acc_rankedTaxids[acc]=L_TaxidsLineage
                            D_Acc_rankedlineage[acc]=L_TaxidsLineageNames
                    except KeyError:
                        print 'acc, 2txID', acc
                TaxaWithSameNames.append([acc, Acc_taxid_temp]) #write the function
                #print L_TaxidsLineage
                #print L_TaxidsLineageNames
        except KeyError:
            L_DelNodes.append(acc)
    for acc_tx in D_Acc_rankedTaxids:
        lastTaxon_taxid=D_Acc_rankedTaxids[acc_tx][-1].split('_')[-1]
        try:
            mCode=Taxid_MitoCode[lastTaxon_taxid]
            D_MitoCode_Acc.setdefault(mCode, []).append(acc_tx)
        except KeyError:
            print '1->', acc_tx
    if len(L_DelNodes)>0:
        f_out_checkDelNodes=open(GoodIds_f.split('.')[0]+'_CheckDelNodes.txt','w')
        f_out_checkDelNodes.write('\n'.join(L_DelNodes)+'\n')
        f_out_checkDelNodes.close()
    for AccNotF in SciNamesNotFound_HigherTaxID: #this is to recognize the MitoCode for accs not found in taxonomy, so we mapped the higher level taxon
        try:
            txId=SciNamesNotFound_HigherTaxID[AccNotF][0]
            MitoCode_HigherTx=Taxid_MitoCode[txId]
            D_MitoCode_Acc.setdefault(MitoCode_HigherTx, []).append(AccNotF)
        except KeyError:
            print '2->', AccNotF       
    return D_Acc_rankedlineage, D_Acc_rankedTaxids, D_MitoCode_Acc, TaxaWithSameNames



SelectedAccs_rankedLineage=get_SelectedAccs_rankedLineage(Id_UniqueTaxID, Taxid_lineageNames, Taxid_Lineage, Taxid_Ranks, TaxRanks, Taxid_MitoCode, SciNamesNotFound_HigherTaxID, BOLD_Acc_AllTaxonomyNoNA)
E=SelectedAccs_rankedLineage[0]
F=SelectedAccs_rankedLineage[1]
print F.items()[:20]
D=SelectedAccs_rankedLineage[2]
TaxaWithSameName=SelectedAccs_rankedLineage[3]
print TaxaWithSameName


def writeAcc_Mitocode(MitoCode_Acc):
    for c in MitoCode_Acc:
        if len(MitoCode_Acc[c])>0:
            f_out=open(GoodIds_f.split('.')[0]+'_'+c+'_withMitocode.txt','w')
            f_out.write('\n'.join(list(set(MitoCode_Acc[c]))))
            f_out.close()
            
Mcode=writeAcc_Mitocode(D)


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


