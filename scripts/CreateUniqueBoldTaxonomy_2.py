import os,sys
from Bio import SeqIO

FinalBOLD_Taxfile='FinalBOLDRecs_taxonomy.tsv'


Original_BOLDTaxonomy='AllBoldRecs_COXI_Taxonomy.tsv'

FinalJoinedTaxonomy='BOLD_ENA_Dereplicated_COXI_Taxonomy.tsv'

def format_NCBI_TaxonomyRecs(FinalNCBI_Taxfile):
    D_Acc_Taxonomy={}
    f=open(FinalNCBI_Taxfile,'r')
    for line in f:
        l=line.strip().split('\t')
        Acc=l[0]
        A_Taxonomy=[]
        for t in l[1:]:
            A_Taxonomy.append(t.split('_')[1])
        D_Acc_Taxonomy[Acc]='\t'.join(A_Taxonomy)
    f.close()
    return D_Acc_Taxonomy


#NCBI_Acc_Taxonomy=format_NCBI_TaxonomyRecs(FinalNCBI_Taxfile)

def format_BOLD_TaxonomyRecs(FinalBOLD_Taxfile):
    D_Acc_Taxonomy={}
    f=open(FinalBOLD_Taxfile,'r')
    for line in f:
        l=line.strip().split('\t')
        Acc=l[0]
        A_Taxonomy=[]
        for t in l[1:]:
            t=t.replace('_',' ')
            A_Taxonomy.append(t)
        D_Acc_Taxonomy[Acc]='\t'.join(A_Taxonomy)
    f.close()
    return D_Acc_Taxonomy

BOLD_Acc_Taxonomy=format_BOLD_TaxonomyRecs(FinalBOLD_Taxfile)

BOLD_Acc_Taxonomy_Original=format_BOLD_TaxonomyRecs(Original_BOLDTaxonomy)

JoinedTaxonomy=format_BOLD_TaxonomyRecs(FinalJoinedTaxonomy)


def get_Phylum2Species(Str_taxonomy): #in the original BOLD we have not the kingdom name, so we remove it in the new generated taxonomy to compare them. 
    taxo_temp='\t'.join(Str_taxonomy.split('\t')[1:])
    return taxo_temp



def getUniqueBoldTaxonomy(BOLD_Acc_Taxonomy_Original, BOLD_Acc_Taxonomy, JoinedTaxonomy):
    D_Acc_Taxonomy_Unique={}
    D_Acc_Taxonomy_Changed={}
    D_Acc_Taxonomy_UnChanged={}
    for acc in BOLD_Acc_Taxonomy_Original:
        if JoinedTaxonomy.has_key(acc): 
            if BOLD_Acc_Taxonomy_Original[acc] == get_Phylum2Species(JoinedTaxonomy[acc]):
                D_Acc_Taxonomy_UnChanged[acc]=JoinedTaxonomy[acc]
            else: #BOLD_Acc_Taxonomy_Original[acc] != get_Phylum2Species(JoinedTaxonomy[acc]):
                D_Acc_Taxonomy_Changed[acc]=[JoinedTaxonomy[acc],BOLD_Acc_Taxonomy_Original[acc]]
#             elif BOLD_Acc_Taxonomy.has_key(acc):
#                 assert get_Phylum2Species(JoinedTaxonomy[acc])==BOLD_Acc_Taxonomy_Original[acc]
#                 D_Acc_Taxonomy_Unique[acc]=JoinedTaxonomy[acc]
        #except KeyError:
            #continue
    f_out_core=open('UnchangedBOLD_Taxonomy.tsv','w')
    for u in D_Acc_Taxonomy_UnChanged:
        f_out_core.write(u+'\t'+D_Acc_Taxonomy_UnChanged[u]+'\n')
    #for c in D_Acc_Taxonomy_Changed:
        #f_out_core.write(c+'\t'+D_Acc_Taxonomy_Changed[c][0]+'\n')
    f_out_core.close()
    f_out_updated_recs=open('UpdatedBOLD_Taxonomy.tsv', 'w')
    f_out_updated_recs.write('Accession\tNCBI_Taxonomy\t\tBOLD_Taxonomy\n')
    for c1 in D_Acc_Taxonomy_Changed:
        f_out_updated_recs.write(c1+'\t'+'\t\t'.join(D_Acc_Taxonomy_Changed[c1])+'\n')
    f_out_updated_recs.close()
            
UU=getUniqueBoldTaxonomy(BOLD_Acc_Taxonomy_Original, BOLD_Acc_Taxonomy, JoinedTaxonomy)


