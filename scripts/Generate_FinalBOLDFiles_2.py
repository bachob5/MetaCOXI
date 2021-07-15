import os,sys
from Bio import SeqIO


CompleteFastaFile=sys.argv[1]
NCBItaxonomyMapped=sys.argv[2]
ResidualBOLDtaxonomy=sys.argv[3]
IDsPassedPFAM=sys.argv[4]
NCBItaxonomyMapped_TaxIDs=sys.argv[5]
ResidualBOLDtaxonomy_TaxIDs=sys.argv[6]

def ID_SEQ(CompleteFastaFile):
    D_Id_Seq={}
    fFas=SeqIO.parse(CompleteFastaFile, 'fasta')
    for rec in fFas:
        D_Id_Seq[rec.id]=str(rec.seq)
    return D_Id_Seq

IDsSeqs=ID_SEQ(CompleteFastaFile)


def residual_BOLDTaxo(ResidualBOLDtaxonomy):
    D_Id_BoldTaxo={}
    f=open(ResidualBOLDtaxonomy, 'r')
    for line in f:
        l=line.strip().split('\t')
        Id=l[0]
        BoldTaxo='\t'.join(l[1:])
        D_Id_BoldTaxo[Id]=BoldTaxo
    f.close()
    return D_Id_BoldTaxo

BoldTaxo_Residual=residual_BOLDTaxo(ResidualBOLDtaxonomy)
BoldTaxIDs_Residual=residual_BOLDTaxo(ResidualBOLDtaxonomy_TaxIDs)

def Id_Taxonomy(f_taxo):
    D_Id_Taxo={}
    f=open(f_taxo, 'r')
    for line in f:
        l=line.strip().split('\t')
        Id=l[0]
        Taxo='\t'.join([x.split('_')[1] for x in l[1:]])
        D_Id_Taxo[Id]=Taxo
    f.close()
    return D_Id_Taxo
    
NCBITaxo=Id_Taxonomy(NCBItaxonomyMapped)
NCBITaxo_TaxIDs=Id_Taxonomy(NCBItaxonomyMapped_TaxIDs)



def get_GoodIds(f_PFAMout):
    GoodIDs_StartEnd={}
    f=open(f_PFAMout, 'r')
    for line in f:
        l=line.split('\t')
        Id=l[0]
        St=int(l[1])
        Ed=int(l[2])
        GoodIDs_StartEnd[Id]=[St, Ed]
    f.close()
    return GoodIDs_StartEnd
    
PassedPFAMIds=get_GoodIds(IDsPassedPFAM)



def JoinAll(IDsSeqs, PassedPFAMIds, BoldTaxo_Residual, NCBITaxo, NCBITaxo_TaxIDs, BoldTaxIDs_Residual):
    f_out_fas=open('FinalBOLDRecs.fasta', 'w')
    f_out_taxo=open('FinalBOLDRecs_taxonomy.tsv', 'w')
    f_out_taxids=open('FinalBOLDRecs_taxonomy_TaxIDs.tsv', 'w')
    for Id in PassedPFAMIds:
        St=PassedPFAMIds[Id][0]
        Ed=PassedPFAMIds[Id][1]
        try:
            f_out_taxo.write(Id+'\t'+NCBITaxo[Id]+'\n')
            f_out_taxids.write(Id+'\t'+NCBITaxo_TaxIDs[Id]+'\n')
            f_out_fas.write('>'+Id+'\n'+IDsSeqs[Id][St:Ed]+'\n')
        except KeyError:
            try:
                f_out_taxo.write(Id+'\t'+BoldTaxo_Residual[Id]+'\n')
                f_out_taxids.write(Id+'\t'+BoldTaxIDs_Residual[Id]+'\n')
                f_out_fas.write('>'+Id+'\n'+IDsSeqs[Id][St:Ed]+'\n')
            except KeyError:
                print Id
    f_out_fas.close()
    f_out_taxo.close()
    f_out_taxids.close()


JoinAll(IDsSeqs, PassedPFAMIds, BoldTaxo_Residual, NCBITaxo, NCBITaxo_TaxIDs, BoldTaxIDs_Residual)





