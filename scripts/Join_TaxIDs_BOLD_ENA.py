import sys

FinalTaxonomyFile=sys.argv[1]
BOLD_TaxIds=sys.argv[2]
ENA_TaxIds=sys.argv[3]

def getFinalIDs(FinalTaxonomyFile):
    L_IDs=[]
    f=open(FinalTaxonomyFile, 'r')
    for line in f:
        l=line.split('\t')
        L_IDs.append(l[0])
    f.close()
    return L_IDs

Ids=getFinalIDs(FinalTaxonomyFile)

def format_ENA_TaxidsRecs(FinalENA_Taxfile):
    D_Acc_TaxIds={}
    f=open(FinalENA_Taxfile,'r')
    for line in f:
        l=line.strip().split('\t')
        Acc=l[0]
        A_Taxids=[]
        for t in l[1:]:
            A_Taxids.append(t.split('_')[1])
        D_Acc_TaxIds[Acc]='\t'.join(A_Taxids)
    return D_Acc_TaxIds

ENA_TaxIds_Recs=format_ENA_TaxidsRecs(ENA_TaxIds)


def format_BOLD_TaxidsRecs(FinalBOLD_Taxfile):
    D_Acc_Taxids={}
    f=open(FinalBOLD_Taxfile,'r')
    for line in f:
        l=line.strip().split('\t')
        Acc=l[0]
        A_Taxids=[]
        for t in l[1:]:
            A_Taxids.append(t)
        D_Acc_Taxids[Acc]='\t'.join(A_Taxids)
    return D_Acc_Taxids

BOLD_Taids_Recs=format_BOLD_TaxidsRecs(BOLD_TaxIds)

f_out=open('BOLD_ENA_Dereplicated_COXI_Taxonomy_Taxids.tsv','w')
for Id in Ids:
    try:
        ENA_TaxIds_Recs[Id]
        f_out.write(Id+'\t'+ENA_TaxIds_Recs[Id]+'\n')
    except KeyError:
        f_out.write(Id+'\t'+BOLD_Taids_Recs[Id]+'\n')

f_out.close()

