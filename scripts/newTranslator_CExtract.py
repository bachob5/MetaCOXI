import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import re
import math
from Bio.Data.CodonTable import TranslationError
"""
MadTranslator.py translate nucleotide sequence in protein for arbitrary genetic codes and frames
and add a suffix to the name of the sequence such that different translation of the same sequence
could be distinguished. The suffix is of the form "seqname_codeX_frameY" where X is one
of the integer recgnized by NCBI for genetic code and Y  is the frame using the standard explained below

python MadTranslator.py inputfasta codesList framelist

inputfasta path of the  multifasta to be translated
codesList  string using the ncbi number genetic code es. "1,2,3"
framelist  string using the notation to identify the six possible frames es. "1,2,3,-1,-2,-3"

restfull interface
inputfasta codesList framelist

example file AraneaeNuc.fas
example command
AraneaeNuc.fas 1,4,5 1,-1,2,3

"""

def startSeqCalc(X, Y, F, LseqOrg, LseqProt):
    """X is the start position. All counts start from '0'; Lseq is the length of the original DNA sequence"""
    if F==1:
        Start = ((X+1)*3)-3
        End = Y*3
        return Start, End
    elif F== -1:
        Start = LseqOrg-(X*3)
        #Start=End- (LseqProt*3)
        End=LseqOrg-(Y*3) -1
        return Start, End
    elif F==2:
        Start = ((X+1)*3)-2
        End = Y*3 +1
        return Start, End
    elif F== -2:
        Start = LseqOrg-(X*3) -1
        End=LseqOrg-(Y*3) -2
        return Start, End
    elif F==3:
        Start = ((X+1)*3)-1
        End = Y*3 +2
        return Start, End
    elif F== -3:
        Start = LseqOrg-(X*3) -2
        End=LseqOrg-(Y*3) -3
        return Start, End

def FlankSLeftFor(F):
    if F ==1:
        return 0
    elif F == 2:
        return 1
    elif F ==3:
        return 2

def FlankSRightFor(X,F):
    if F ==1:
        return ((X+1)*3)
    elif F == 2:
        return ((X+1)*3) +1
    elif F ==3:
        return ((X+1)*3) +2

def FlankSLeftRev(X, F, LseqOrg):
    if F ==-1:
        Start = LseqOrg - (X*3)
        End = LseqOrg
        return Start, End
    elif F == -2:
        Start = LseqOrg - (X*3) - 1
        End = LseqOrg -1
        return Start, End
    elif F == -3:
        Start = LseqOrg - (X*3) - 2
        End = LseqOrg -2
        return Start, End

def FlankSRightRev(X, F, LseqOrg):
    if F ==-1:
        End = LseqOrg - (X*3)
        return End
    elif F == -2:
        End = LseqOrg - (X*3) -1
        return End
    elif F == -3:
        End = LseqOrg - (X*3) -2
        return End

#def startEndSeqCalc_Complement(X, F, LseqOrg, LseqProt):

def getFrag_StopCodons(Sequence, ID, frm, LenSeq, MinValue=5):
    Dict={}
    A=re.finditer('@', Sequence)
    B=[x.span() for x in A]
    T=[]
    for n in B:
        T.append(list(n))
    #print T
    C=sum(T,[])
    #print C
    if len(C)>0:
        if int(frm) > 0:
            if C[0] > MinValue:
                Start = FlankSLeftFor(int(frm))
                newID=ID+';S'+ str(Start)+ ';E'+ str((C[0]*3)+(int(frm)-1))
                #print newID, C[0]
                Dict[newID] = Sequence[:C[0]]
            if len(Sequence)-C[-1] > MinValue:
                Start=FlankSRightFor(C[-1], int(frm))
                newID=ID+';S' + str(Start)+';E'+ str((len(Sequence)*3)+int(frm)-1)
                #print newID, C[-1]
                Dict[newID] = Sequence[C[-1]+1:]
            i=0
            while i<len(C)-3:
                if C[i+2]-C[i+1] > MinValue:
                    SliceProt=C[i+2]-C[i+1]
                    start=startSeqCalc(C[i+1], C[i+2], int(frm), LenSeq, SliceProt)
                    newID=ID+';S' + str(start[0]) + ';E'+ str(start[1])
                    #print newID, C[i+1], C[i+2]
                    Dict[newID] = Sequence[C[i+1]:C[i+2]]
                i=i+2
        else:
            if C[0] > MinValue:
                StartEn = FlankSLeftRev(C[0], int(frm), LenSeq)
                newID=ID+';S' + str(StartEn[0])+';E' + str(StartEn[1])
                Dict[newID] = Sequence[:C[0]]
                #print newID, C[0]
            if len(Sequence)-C[-1] > MinValue:
                #St=LenSeq-(C[-1]*3)
                St='0'
                En=FlankSRightRev(C[-1], int(frm), LenSeq)
                strt= En-((len(Sequence) - C[-1])*3)
                newID=ID+';S' + str(strt)+';E'+ str(En)
                #print newID, C[-1], LenSeq, len(Sequence)
                Dict[newID] = Sequence[C[-1]:]
            j=0
            while j<len(C)-3:
                if C[j+2]-C[j+1] > MinValue:
                    SliceProt=C[j+2]-C[j+1]
                    startEnd=startSeqCalc(C[j+2], C[j+1], int(frm), LenSeq, SliceProt)
                    #print startEnd[0]
                    newID=ID+';S' + str(startEnd[0]) + ';E'+ str(startEnd[1]+1)
                    #print newID, C[j+1], C[j+2]
                    Dict[newID] = Sequence[C[j+1]:(C[j+2])]
                j=j+2
    elif len(C) <= 0: 
        if int(frm) > 0:
            if frm == '1':
                newID=ID+';S0' + ';E'+ str((len(Sequence)*3)+1)
                #print newID
                Dict[newID] = Sequence
            elif frm =='2':
                newID=ID+';S1' + ';E'+ str((len(Sequence)*3)+1)
                #print newID
                Dict[newID] = Sequence
            elif frm =='3':
                newID=ID+';S2' + ';E'+ str((len(Sequence)*3)+1)
                #print newID
                Dict[newID] = Sequence
        else:
            St='0'
            if frm == '-1':
                End=LenSeq
                Strt=LenSeq-len(Sequence)*3
                newID=ID+';S' + str(Strt)+';E'+ str(End)
                #print newID, C
                Dict[newID] = Sequence
            elif frm == '-2':
                End=LenSeq-1
                Strt=LenSeq-len(Sequence)*3
                newID=ID+';S' + str(Strt)+';E'+ str(End)
                #print newID, C
                Dict[newID] = Sequence
            elif frm == '-3':
                End=LenSeq-2
                Strt=LenSeq-len(Sequence)*3
                newID=ID+';S' + str(Strt)+';E'+ str(End)
                #print newID, C
                Dict[newID] = Sequence
    return Dict

codes=[5,2,3,4,9,10,11,12,13,14,15,1]
codes=map(int,sys.argv[2].split(','))
tempframe=sys.argv[3].split(',')
frames=zip([abs(int(x))-1 for x in tempframe ],[ int(x)<0 for x in tempframe],tempframe)
#print frames
fName=sys.argv[1]
SeqGen=SeqIO.parse(fName,'fasta')
seqs={}
f_out=open(fName.split('.')[0]+'_codeframes.fas', 'w')
for seq in SeqGen:
    name=seq.id
    LenSeq1=len(seq.seq)
    for t in codes:
        for frame in frames:
            try:
                if frame[1]:
                    Seq=seq.reverse_complement()
                    temp=str(Seq[frame[0]:].seq.translate(table=t, stop_symbol="@"))
                else:
                    temp=str(seq[frame[0]:].seq.translate(table=t, stop_symbol="@"))
                newName = name+';code'+str(t)+';frame'+str(frame[2])
                seqs=getFrag_StopCodons(temp, newName, frame[2], LenSeq1)
                for ID in seqs:
                    f_out.write('>' + ID + '\n' + seqs[ID]+'\n')
            except TranslationError:
                continue
            #seqs.append(SeqRecord(seq=temp,description='',id=name+'_code'+str(t)+'_frame'+str(frame[2])))

f_out.close()
