import os, sys
import re
from Bio import SeqIO

f_fas=sys.argv[1]
SequenceLen=sys.argv[2]
NumberNs=sys.argv[3]

f_delnodes=sys.argv[4]
f_notFoundinTax=sys.argv[5]
f_PseudoIds=sys.argv[6]

def getAllLines(f_in):
    Ids=[]
    if os.path.isfile(f_in):
        f=open(f_in, 'r')
        for line in f:
            Ids.append(line.strip())
    else:
        pass
    return Ids


DelNodes=getAllLines(f_delnodes)
NotInTax=getAllLines(f_notFoundinTax)
PseudoIds=getAllLines(f_PseudoIds)

IdsToDel=DelNodes+NotInTax+PseudoIds
#print IdsToDel

def transformList2Dict(L):
    D_Ids={}
    for i in L:
        D_Ids[i]=True
    return D_Ids

D_Ids2Del=transformList2Dict(IdsToDel)

def get_IdsFromTaxFile(f_fas, Ext, D_Ids2Del):
    D_IdsTaxonomy={}
    f_in=open(f_fas.split('.')[0]+'_'+Ext, 'r')
    for line in f_in:
        l=line.strip().split('\t')
        try:
            D_Ids2Del[l[0]]
        except KeyError:
            D_IdsTaxonomy[l[0]]=True
    f_in.close()
    return D_IdsTaxonomy

D_IdsTaxo=get_IdsFromTaxFile(f_fas, 'Taxonomy_TaxIDs.txt', D_Ids2Del)


def trim_Ns_Ends(DNAseq):
    S={}
    E={}
    Ns2Trim_S_E=[]
    NNs=re.findall('N+', DNAseq)
    for n in NNs:
        if DNAseq.startswith(n):
            S.setdefault(n,[]).append(len(n))
        if DNAseq.endswith(n):
            E.setdefault(n,[]).append(len(n))
    #for pat in S:
    #    if len(S[pat])==2:
    #        assert E[pat]==2
    #        raise AssertionError
    if len(S)>0:
        LensS=sum(S.values(),[])
        LensS.sort(reverse=True)
        Ns2Trim_S_E.append(LensS[0])
    else:
        Ns2Trim_S_E.append(0)
    if len(E)>0:
        LensE=sum(E.values(),[])
        LensE.sort(reverse=True)
        Ns2Trim_S_E.append(len(DNAseq) - LensE[0])
    else:
        Ns2Trim_S_E.append(len(DNAseq))
    return Ns2Trim_S_E



def writeFinalSeqFileWithoutNNs(f_fas, D_IdsTaxo):
    DremoveIds={}
    D_IdsToKeep={}
    D_Id_newDescription={}
    fas=SeqIO.parse(f_fas, 'fasta')
    f_out=open(f_fas.split('.')[0]+'_WithoutNs.fasta', 'w')
    for rec in fas:
        try:
            D_IdsTaxo[rec.id]
            DNAseq=str(rec.seq)
            Ends2Trim=trim_Ns_Ends(DNAseq)
            S=Ends2Trim[0]
            E=Ends2Trim[1]
            newSeq=DNAseq[S:E]
            oldDescription=rec.description.split()
            oldPoss=map(int, oldDescription[1].split('-'))
            if newSeq.count('N')<=int(NumberNs) and len(newSeq)>=int(SequenceLen):
                #if oldPoss==Ends2Trim:
                #    f_out.write('>'+rec.description+'\n'+DNAseq+'\n')
                #else:
                newPoss=oldDescription[1]+'@'+str(S)+'-'+str(E)
                newDescription=rec.id+' '+newPoss+' '+oldDescription[2]
                f_out.write('>'+newDescription+'\n'+newSeq+'\n')
            else:
                DremoveIds[rec.id]=True
        except KeyError:
            continue
    f_out.close()
    #LremoveIds=list(set(LremoveIds_temp+IdsToDel))
    for Id in D_IdsTaxo:
        try:
            DremoveIds[Id]
        except KeyError:
            D_IdsToKeep[Id]=True
    return D_IdsToKeep

D_Ids2Keep=writeFinalSeqFileWithoutNNs(f_fas, D_IdsTaxo)
print len(D_Ids2Keep)


def createNewTaxFile(f_fas, D_Ids2Keep, Ext):
    fName_base=f_fas.split('.')[0]+'_'+Ext
    f_in=open(f_fas.split('.')[0]+'_'+Ext, 'r')
    f_out=open(fName_base.split('.')[0]+'_WithoutNNs.tsv', 'w')
    for line in f_in:
        l=line.strip().split('\t')
        try:
            D_Ids2Keep[l[0]]
            f_out.write(line.strip()+'\n')
        except KeyError:
            continue
    f_in.close()
    f_out.close()

createNewTaxFile(f_fas, D_Ids2Keep, 'Taxonomy_Lineages.txt')
createNewTaxFile(f_fas, D_Ids2Keep, 'Taxonomy_TaxIDs.txt')


