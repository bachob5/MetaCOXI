import sys,os

f_OnlyAccs=sys.argv[1]
f_BoldTaxonomyFile=sys.argv[2]


def get_Accs2Analyze(f_OnlyAccs):
    Accs=[]
    l=[]
    f_in=open(f_OnlyAccs,'r')
    for line in f_in:
        if len(line.strip())!=0:
            Accs.append(line.strip())
        else:
            l.append('line')
    f_in.close()
    print l
    return Accs

Accs2Analyze=get_Accs2Analyze(f_OnlyAccs)
print len(Accs2Analyze)


def get_BOLDTaxNames(f_BoldTaxonomyFile, Accs2Analyze):
    D_Acc_Taxonomy={}
    A_Accs2Analyze_Taxonomy={}
    #L=[]
    f=open(f_BoldTaxonomyFile, 'r')
    for line in f:
        l=line.strip().split('\t')
        D_Acc_Taxonomy[l[0]]=line.strip()
        #D_Acc_Taxonomy.setdefault(l[0],[]).append(line.strip())
        #L.append(l[0])
    f.close()
#     for i in D_Acc_Taxonomy:
#         if len(D_Acc_Taxonomy[i])>1:
#             print D_Acc_Taxonomy[i]
#     print len(D_Acc_Taxonomy)
#     print len(L), len(set(L))
    for acc in Accs2Analyze:
        A_Accs2Analyze_Taxonomy[acc]=D_Acc_Taxonomy[acc]
    f_out=open(f_OnlyAccs.split('.')[0]+'_Accs2AnalyzeTaxonomy.txt', 'w')
    for Acc in A_Accs2Analyze_Taxonomy:
        f_out.write(A_Accs2Analyze_Taxonomy[Acc]+'\n')
    f_out.close()

get_BOLDTaxNames(f_BoldTaxonomyFile, Accs2Analyze)