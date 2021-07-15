import sys,os
from Bio import SeqIO

CurrFolder=os.getcwd().split('/')[-1]


PfamSingleMatchIDsFiles=[x for x in os.listdir(os.curdir) if '_'.join(x.split('_')[6:]) =='withMitocode_codeframes_MatchIDs_PfamProfile.txt']
print PfamSingleMatchIDsFiles



def get_PfamAnnotations(f_IDs):
    D_Acc_Cox1Poss={}
    f=open(f_IDs, 'r')
    for line in f:
        l=line.strip().split()
        for i in l:
            temp=i.split(';')
            Id=temp[0]
            Start_Pos=int(temp[-2])
            End_Pos=int(temp[-1])
            D_Acc_Cox1Poss.setdefault(Id,[]).append([Start_Pos, End_Pos])
    f.close()
    return D_Acc_Cox1Poss

# D_PfamAnnotation=get_PfamAnnotations(try_PfamMatchFiles)
# print len(D_PfamAnnotation)

def get_newPoss(oldPos):
    if oldPos.startswith('<'):
        newPos=oldPos[1:]
        return newPos
    elif oldPos.startswith('>'):
        newPos=oldPos[1:]
        return newPos
    else:
        newPos=oldPos
        return newPos
    

def get_GbAnnotations(f_gb):
    D_Acc_Annotations={}
    handle_GB=SeqIO.parse(f_gb, 'genbank')
    for rec in handle_GB:
        tempL=[]
        for seq_feature in rec.features:
            if seq_feature.type=='CDS': # Add the control of Pseudogenes
                if len(seq_feature.location.parts) > 1:
                    tempJoin=[]
                    for j in seq_feature.location.parts:
                        if 'gene' in seq_feature.qualifiers:
                            GeneName=seq_feature.qualifiers['gene'][0]
                            startPos=get_newPoss(str(j.start))
                            endPos=get_newPoss(str(j.end))
                            tempJoin.append('Gene='+GeneName+';'+startPos+';'+endPos)
                        else:
                            if 'product' in seq_feature.qualifiers:
                                ProductName=seq_feature.qualifiers['product'][0]
                                startPos=get_newPoss(str(j.start))
                                endPos=get_newPoss(str(j.end))
                                tempJoin.append('product='+ProductName+';'+startPos+';'+endPos) #this means that there is no annotation named 'gene' in the entry
                    tempL.append(tempJoin)
                elif len(seq_feature.location.parts) == 1: # Add the control of Pseudogenes
                    PosS = str(seq_feature.location.start)
                    newPosS=get_newPoss(PosS)
                    PosE=str(seq_feature.location.end)
                    newPosE=get_newPoss(PosE)
                    if 'gene' in seq_feature.qualifiers:
                        GeneName=seq_feature.qualifiers['gene'][0]
                        tempL.append(['Gene='+GeneName+';'+newPosS+';'+newPosE])
                    else:
                        if 'product' in seq_feature.qualifiers:
                            ProductName=seq_feature.qualifiers['product'][0]
                            tempL.append(['product='+ProductName+';'+newPosS+';'+newPosE])
                                        #print rec.id
        D_Acc_Annotations.setdefault(rec.id, []).append(tempL)
    return D_Acc_Annotations



def get_Intersections_LL(L1, L2):
    range_L1=range(L1[0], L1[1])
    range_L2=range(L2[0], L2[1])
    L_Intersect=set(range_L1).intersection(set(range_L2))
    return list(L_Intersect)


def get_PositiveCOX1(D_PfamAnnotation, D_GbAnnotations):
    D_Id_Gb_Poss_GbAnnotat={}
    L_MissedAccs=[]
    for acc in D_GbAnnotations:
        GB_Rec=D_GbAnnotations[acc]
        for ann in GB_Rec:
            for j in ann:
                for i in j:
                    try:
                        start_annot=int(i.split(';')[-2])
                        end_annot=int(i.split(';')[-1])
                        p_gb=[start_annot, end_annot]
                        try:
                            for Pos in D_PfamAnnotation[acc]:
                                if len(get_Intersections_LL(p_gb, Pos))>30:
                                    D_Id_Gb_Poss_GbAnnotat.setdefault(acc, []).append([p_gb,Pos])
                        except KeyError:
                            continue
                    except IndexError:
                        continue
        if acc not in D_PfamAnnotation:
            L_MissedAccs.append(acc)
        if acc not in D_Id_Gb_Poss_GbAnnotat:
            L_MissedAccs.append(acc)
    return D_Id_Gb_Poss_GbAnnotat, list(set(L_MissedAccs))


def get_COXI_synonyms(D_GbAnnotations, D_Id_Gb_Poss_GbAnnotat):
    L_acc_MultipleMatch=[]
    L_acc_NoMatch=[]
    D_GB_Annotations_GbPoss_PfamPoss={}
    for acc in D_Id_Gb_Poss_GbAnnotat:
        try:
            Poss_gb_pfam=D_Id_Gb_Poss_GbAnnotat[acc]
            if len(Poss_gb_pfam) > 1:
                L_acc_MultipleMatch.append(acc)
            elif len(Poss_gb_pfam)==0:
                L_acc_NoMatch.append(acc)
            elif len(Poss_gb_pfam)==1:
                temp_Poss=map(str,Poss_gb_pfam[0][0])
                temp_Poss_Pfam=map(str, Poss_gb_pfam[0][1])
                Org_annotations=D_GbAnnotations[acc][0]
                for anns in Org_annotations:
                    for ann in anns:
                        temp_Id=ann.split(';')
                        if temp_Id[-2:]==temp_Poss:
                            D_GB_Annotations_GbPoss_PfamPoss.setdefault(acc, []).append([temp_Id[0], temp_Poss, temp_Poss_Pfam])
        except IndexError:
            continue
    return D_GB_Annotations_GbPoss_PfamPoss, L_acc_MultipleMatch, L_acc_NoMatch



def write_GoodCOX1Seqs(D_GB_Annotations_GbPoss_PfamPoss, f_gb):
    D_Acc_Poss={}
    L_COX1_synonyms={}
    L_Acc_withPotentialFrameShift=[]
    f_out=open(f_gb.split('.')[0]+'_SelectedGoodCOX1.fasta', 'w')
    for acc in D_GB_Annotations_GbPoss_PfamPoss:
        for p in D_GB_Annotations_GbPoss_PfamPoss[acc]:
            L_COX1_synonyms.setdefault(acc,[]).append(p[0])
            gb_poss=map(int, p[1])
            D_Acc_Poss.setdefault(acc,[]).append(gb_poss)
    f=SeqIO.parse(f_gb, 'genbank')
    for rec in f:
        try:
            if len(D_Acc_Poss[rec.id])>1:
                L_Acc_withPotentialFrameShift.append(rec.id) #write the function later
            elif len(D_Acc_Poss[rec.id])==1:
                start_pos=D_Acc_Poss[rec.id][0][0]
                end_pos=D_Acc_Poss[rec.id][0][1]
                DNAseq=str(rec.seq)[start_pos:end_pos]
                GeneName=''.join(L_COX1_synonyms[rec.id])
                if len(DNAseq)>100:
                    f_out.write('>'+rec.id+' '+str(start_pos)+'-'+str(end_pos)+' '+ GeneName +'\n'+DNAseq+'\n')
        except KeyError:
            continue
    f_out.close()
    temp_COXISynonyms=[] #list(set(L_COX1_synonyms))
    for a in L_COX1_synonyms:
        temp_COXISynonyms.append(L_COX1_synonyms[a])
    COXISynonyms=list(set(sum(temp_COXISynonyms,[])))
    if len(COXISynonyms)>0:
        f_out_COXI_Synonyms=open(f_gb.split('.')[0]+'_COXISynonyms.txt', 'w')
        f_out_COXI_Synonyms.write('\n'.join(COXISynonyms)+'\n')
        f_out_COXI_Synonyms.close()
    AccsWithFrameShifts=list(set(L_Acc_withPotentialFrameShift))
    if len(AccsWithFrameShifts)>0:
        f_out_FrameShifts=open(f_gb.split('.')[0]+'_AccsWithPotentialFrameShifts.txt', 'w')
        f_out_FrameShifts.write('\n'.join(AccsWithFrameShifts)+'\n')
        f_out_FrameShifts.close()


#Res_GoodCOXI=write_GoodCOX1Seqs(D_GB_Annotations_GbPoss_PfamPoss, try_GBFiles)


def get_MissedAccessions(MissedAccs, L_acc_MultipleMatch, L_acc_NoMatch, f_gb):
    MissedAll=MissedAccs+L_acc_MultipleMatch+L_acc_NoMatch
    if len(MissedAll)>0:
        f_out=open(f_gb.split('.')[0]+'_VerifyMissedAccs.txt', 'w')
        f_out.write('\n'.join(MissedAll))
        f_out.close()


#M=get_MissedAccessions(MissedAccs, L_acc_MultipleMatch, L_acc_NoMatch, try_GBFiles)


for PfamMatchFile in PfamSingleMatchIDsFiles:
    pf_f_temp=PfamMatchFile.split('_')
    GBFile='_'.join(pf_f_temp[:5])+'_Good_Excluded_Entries.gb' #remove '_SelectedEntries_' in case u use raw gb files
    D_PfamAnnotation=get_PfamAnnotations(PfamMatchFile)
    print PfamMatchFile, len(D_PfamAnnotation)
    D_GbAnnotations=get_GbAnnotations(GBFile)
    print GBFile, len(D_GbAnnotations)
    Annotations = get_PositiveCOX1(D_PfamAnnotation, D_GbAnnotations)
    D_Id_Gb_Poss_GbAnnotat=Annotations[0]
    #print D_Id_Gb_Poss_GbAnnotat
    MissedAccs=Annotations[1]
    print 'D_Id_Gb_Poss_GbAnnotat', len(D_Id_Gb_Poss_GbAnnotat) 
    Res_All=get_COXI_synonyms(D_GbAnnotations, D_Id_Gb_Poss_GbAnnotat)
    D_GB_Annotations_GbPoss_PfamPoss=Res_All[0]
    L_acc_MultipleMatch=Res_All[1]
    L_acc_NoMatch=Res_All[2]
    Res_GoodCOXI=write_GoodCOX1Seqs(D_GB_Annotations_GbPoss_PfamPoss, GBFile)
    M=get_MissedAccessions(MissedAccs, L_acc_MultipleMatch, L_acc_NoMatch, GBFile)


