import os,sys
import re

Markers= 'COI-5P,COI-3P'
list_Markers=Markers.split(',')
#print list_Markers

SequenceLen=sys.argv[1]
NumberNs=sys.argv[2]

def isEmpty(S):
    if S.strip()=='':
        return 'NA'
    else:
        return S

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
        

def writeFasta(f_tsv, GenMarker): #GenMarker is a list of markers. The first element in the ID The Taxonomy is in this order [phylum, class, order, family, genus,species]
    f=open(f_tsv, 'r')
    D_ID_Seq={}
    taxName=f_tsv.split('_FullRecords_BOLD.tsv')[0]
    headers=f.readline().strip().split('\t')
    #SeqID_Num=headers.index('sequenceID')
    ProcessID=headers.index('processid')
    NCBIAcc=headers.index('genbank_accession')
    marker_Num=headers.index('markercode')
    Nuc_Num=headers.index('nucleotides')#.strip().upper()
    #Nuc_Num=Nuc_Num.replace('-','')
    phylum_name=headers.index('phylum_name')
    class_name=headers.index('class_name')
    order_name=headers.index('order_name')
    family_name=headers.index('family_name')
    genus_name=headers.index('genus_name')
    species_name=headers.index('species_name')
    f_out=open(GenMarker + '_'+taxName+'_extractedRecsByMarker.fasta', 'w')
    for line in f:
        l=line.split('\t')
        try:
            DNAseq_temp=l[Nuc_Num].replace('-','')
            PotentialTrimPoss=trim_Ns_Ends(DNAseq_temp)
            #print PotentialTrimPoss
            S=PotentialTrimPoss[0]
            E=PotentialTrimPoss[1]
            #print type(S), type(E)
            DNAseq_new=DNAseq_temp[S:E]
            if l[marker_Num] in GenMarker and DNAseq_new.count('N')<=int(NumberNs) and len(DNAseq_new)>= int(SequenceLen):
                if l[NCBIAcc].strip() == '':
                    D_ID_Seq.setdefault(l[ProcessID]+ '.'+ l[marker_Num] ,[]).append(';' +';'.join([isEmpty(l[phylum_name]),isEmpty(l[class_name]), isEmpty(l[order_name]), isEmpty(l[family_name]), isEmpty(l[genus_name]), isEmpty(l[species_name]).replace(' ','_')]) +'\n' + DNAseq_new.replace('-','') + '\n')
                    #f_out.write('>' + l[ProcessID]+ '.'+ l[marker_Num] + ';' +';'.join([isEmpty(l[phylum_name]),isEmpty(l[class_name]), isEmpty(l[order_name]), isEmpty(l[family_name]), isEmpty(l[genus_name]), isEmpty(l[species_name]).replace(' ','_')]) +'\n' + l[Nuc_Num].replace('-','') + '\n')
                elif l[NCBIAcc].startswith('Pending'):
                    D_ID_Seq.setdefault(l[ProcessID]+ '.'+ l[marker_Num] ,[]).append(';' +';'.join([isEmpty(l[phylum_name]),isEmpty(l[class_name]), isEmpty(l[order_name]), isEmpty(l[family_name]), isEmpty(l[genus_name]), isEmpty(l[species_name]).replace(' ','_')]) +'\n' + DNAseq_new.replace('-','') + '\n')
                    #f_out.write('>' + l[ProcessID]+ '.'+ l[marker_Num] + ';' +';'.join([isEmpty(l[phylum_name]),isEmpty(l[class_name]), isEmpty(l[order_name]), isEmpty(l[family_name]), isEmpty(l[genus_name]), isEmpty(l[species_name]).replace(' ','_')]) +'\n' + l[Nuc_Num].replace('-','') + '\n')
                elif 'SUPPRESSED' in l[NCBIAcc]:
                    continue
                else:
                    D_ID_Seq.setdefault(l[NCBIAcc].split()[0] ,[]).append(';'+ l[marker_Num] + ';' +';'.join([isEmpty(l[phylum_name]),isEmpty(l[class_name]), isEmpty(l[order_name]), isEmpty(l[family_name]), isEmpty(l[genus_name]), isEmpty(l[species_name]).replace(' ','_')]) +'\n' + DNAseq_new.replace('-','') + '\n')
                    #f_out.write('>' + l[NCBIAcc].split()[0] + ';'+ l[marker_Num] + ';' +';'.join([isEmpty(l[phylum_name]),isEmpty(l[class_name]), isEmpty(l[order_name]), isEmpty(l[family_name]), isEmpty(l[genus_name]), isEmpty(l[species_name]).replace(' ','_')]) +'\n' + l[Nuc_Num].replace('-','') + '\n')
        except IndexError:
            print l[0]
    for rec in D_ID_Seq:
        if len(set(D_ID_Seq[rec]))==1:
            f_out.write('>'+rec+D_ID_Seq[rec][0])
        else:
            print 'Unresolved\n', f_tsv, rec
    f.close()
    f_out.close()

filesList=[x for x in os.listdir(os.curdir) if '_FullRecords_BOLD.tsv' in x]

for marker in list_Markers:
    for f_in in filesList:
        print 'processing ', f_in
        if os.path.getsize(f_in)>0:
            writeFasta(f_in, marker)

