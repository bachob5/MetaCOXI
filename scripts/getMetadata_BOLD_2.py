import dateutil
from dateutil import parser
import datetime
from datetime import datetime 
import sys


L=['collectiondate_start', 'collectiondate_end', 'collectiontime', 'associated_taxa', 'lat', 'lon', 'country', 'province_state', 'region', 'sector', 'exactsite', 'markercode', 'genbank_accession']
s='processid'

IDs2Search_BOLD=sys.argv[1]
FullBold=sys.argv[2]

Markers=['COI-3P', 'COI-5P']
BOLD_RetrieveLink='http://www.boldsystems.org/index.php/Public_RecordView?processid='

def get_IDs(IDs2Search_BOLD):
    D_IDs={}
    f=open(IDs2Search_BOLD,'r')
    for line in f:
        Id_Bold=line.strip()
        D_IDs[Id_Bold]=True
    f.close()
    return D_IDs

IdsBOLD=get_IDs(IDs2Search_BOLD)

def get_datetime(Str_datetime):
    d_temp=datetime.utcnow()
    YearNow=d_temp.year
    Monthnow=d_temp.strftime('%b')
    Today=d_temp.day
    TodaysDate=str(Today)+'-'+Monthnow+'-'+str(YearNow)
    try:
        x1=dateutil.parser.parse(Str_datetime, yearfirst=True)
        Year=x1.year
        month_abre=x1.strftime('%b')
        Day=x1.day
        DATE=str(Day)+'-'+month_abre+'-'+str(Year)
        if DATE==TodaysDate:
            return 'NA\tNA'
        else:
            return str(Year) + '\t' + DATE
    except ValueError:
        return 'NA\tNA'
    
'''
datetime.datetime(206, 6, 12, 0, 0)
'''
    
def getLength(Str_value):
    if len(Str_value.strip())>0:
        return Str_value.strip()
    else:
        return 'NA'
    
def getHostName(Str_Host):
    Host_temp=Str_Host.split(':')
    if len(Host_temp)>1:
        return Host_temp[1]
    else:
        return getLength(Str_Host)

def getMetadata_BoldFull(FullBold):
    D_variable_index={}
    D_processID_Metadata={}
    D_GenbankAcc_ProcessID={}
    Indexes=[]
    f=open(FullBold,'r')
    headers=f.readline().strip().split('\t')
    print headers
    D_Index_header={}
    for index,name in enumerate(headers):
        if name in L:
            D_variable_index[str(index)]=name
            Indexes.append(index)
    print Indexes
    for line in f:
        l=line.split('\t')
        #L1=[]
        #for ind in Indexes:
            #if ind in Indexes:
        #collectionDate=l[Indexes[0]]+' '+l[Indexes[1]]+' '+l[Indexes[2]]
        collectionDate_temp=l[Indexes[2]]
        collectionDate=get_datetime(collectionDate_temp)
        host_sp_temp=l[Indexes[3]]
        host_sp=getHostName(host_sp_temp)
        GeoCoor_temp=l[Indexes[4]]+' '+l[Indexes[5]]
        GeoCoor=getLength(GeoCoor_temp)
        country_temp=l[Indexes[6]]
        country=getLength(country_temp)
        c_Address_temp=l[Indexes[7]]+' '+l[Indexes[8]]+' '+l[Indexes[9]]+' '+l[Indexes[10]]
        c_Address=getLength(c_Address_temp)
        marker=l[Indexes[11]]
        #marker=l[Indexes[-2]]
        genbank_acc=getLength(l[Indexes[12]])
        if marker in Markers:
            L1='\t'.join([country, c_Address, host_sp, collectionDate, GeoCoor, genbank_acc])
            D_processID_Metadata.setdefault(l[0]+'.'+marker, []).append(L1)
            if len(genbank_acc)>2:
                D_GenbankAcc_ProcessID[genbank_acc]=l[0]+'.'+marker
    f.close()
    return D_processID_Metadata, D_GenbankAcc_ProcessID


D=getMetadata_BoldFull(FullBold)
ProcessID_Metadata=D[0]
Acc_ProcessID=D[1]

def get_NeededMetadata(IdsBOLD, ProcessID_Metadata, Acc_ProcessID, BOLD_RetrieveLink):
    D_IdsBOLD_Metadata={}
    L_MissedIds=[]
    for Id in IdsBOLD:
        if ProcessID_Metadata.has_key(Id):
            Id_Meta_temp='\t'.join(ProcessID_Metadata[Id][0].split('\t')[:-1])
            Id_Meta=Id_Meta_temp+'\t'+BOLD_RetrieveLink+Id.split('.')[0]
            D_IdsBOLD_Metadata[Id]=Id_Meta
        elif Acc_ProcessID.has_key(Id):
            ProID=Acc_ProcessID[Id]
            try:
                Id_Meta_temp='\t'.join(ProcessID_Metadata[ProID][0].split('\t')[:-1])
                Id_Meta=Id_Meta_temp+'\t'+BOLD_RetrieveLink+ProID.split('.')[0]
                D_IdsBOLD_Metadata[ProID+';'+Id]=Id_Meta
            except KeyError:
                print Id
        else:
            L_MissedIds.append(Id)
    return D_IdsBOLD_Metadata, L_MissedIds

M=get_NeededMetadata(IdsBOLD, ProcessID_Metadata, Acc_ProcessID, BOLD_RetrieveLink)
M_dataBold=M[0]
MissedAccs=M[1]

if len(MissedAccs)>0:
    f_out_MissedAccs=open('MissedAccsVerify.txt', 'w')
    f_out_MissedAccs.write('\n'.join(MissedAccs)+'\n')
    f_out_MissedAccs.close()

f_out=open('Metadata_BOLD.tsv', 'w')
f_out.write('\t'.join(['Accession', 'country', 'country_address' ,'host', 'collection_year' ,'collection_complete_date', 'geographical_coordinates', 'RetrievalLink'])+'\n')
for Acc in M_dataBold:
    f_out.write(Acc+'\t'+M_dataBold[Acc]+'\n')
f_out.close()



'''
['', 'late austral winter', '2015-06-25T18:00:00.000Z', '2016-05-24T17:00:00.000Z', '2017-05-17T13:02:00.000Z', '23:45-4:00', '2017-06-14T12:45:00.000Z', '2017-09-14T23:12:00.000Z', '2016-07-12T17:00:00.000Z', '2016-08-03T13:16:00.000Z', '2016-09-13T16:00:00.000Z', '2017-06-16T01:31:00.000Z', '2016-08-11T11:50:00.000Z', '2017-06-07T11:57:00.000Z', '2017-06-01T20:00:00.000Z', '2016-09-13T15:21:00.000Z', '2017-07-27T10:42:00.000Z', '2017-10-04T13:40:00.000Z', '2017-06-28T16:00:00.000Z', '2017-08-25T13:16:00.000Z', '2016-05-12T16:00:00.000Z', '11:05:00', '2016-07-20T20:00:00.000Z', '07:45-08:00', 'month', '09:30-10:30', '10 AM - 2:00 PM', '2017-08-09T12:43:00.000Z', '2017-08-16T16:00:00.000Z', '2016-05-25T17:00:00.000Z', '2017-11-01T21:52:00.000Z', '2017-07-25T15:24:00.000Z', 'Sunrise', '9pm', '2017-08-24T14:07:00.000Z', 'Daytime', '20:15-22:10', '2016-08-30T22:01:00.000Z', '22:45', '2017-09-21T21:55:00.000Z', '2016-09-08T13:55:00.000Z', '2017-06-20T20:43:00.000Z', '2017-05-16T12:53:00.000Z', '2017-07-20T23:07:00.000Z', '22:48', '15-04-2016', '2017-05-10T12:14:00.000Z', '15:45-16:45', '09:45 -12:00', '2016-10-05T00:00:00.000Z', '2017-04-28T16:00:00.000Z', '2017-05-03T13:04:00.000Z', '14-09-21', '2016-07-29T12:19:00.000Z', '2016-05-09T00:00:00.000Z', '2017-08-15T13:24:00.000Z', '2016-05-10T18:00:00.000Z', '2016-06-17T10:47:00.000Z', '2017-08-30T13:20:00.000Z', '2015-07-29T17:00:00.000Z', '2017-08-28T18:00:00.000Z', '2017-08-24T13:13:00.000Z', '2017-07-18T13:36:00.000Z', '2017-07-18T19:00:00.000Z', '2017-09-26T12:56:00.000Z', '2017-07-14T13:35:00.000Z', '2016-06-30T12:28:00.000Z', '2017-06-20T22:33:00.000Z', '2014-06-11T16:00:00.000Z', '2017-08-09T13:19:00.000Z', '12:19', '2014-06-05T18:00:00.000Z', '2016-05-31T18:00:00.000Z', '19:30-1:30', '2017-08-02T14:11:00.000Z', '2016-09-23T13:00:00.000Z', '2016-06-30T00:00:00.000Z', '2017-09-20T16:00:00.000Z', '2017-08-16T12:56:00.000Z', '2017-07-26T11:29:00.000Z', '8:51', '10:30-11:30', '2017-05-17T12:55:00.000Z', '2016-07-10T00:00:00.000Z', '2014-07-16T16:00:00.000Z', '2016-06-24T13:15:00.000Z', '2017-06-15T11:24:00.000Z', '2016-09-27T15:16:00.000Z', '20:59', '2017-09-20T18:00:00.000Z', '2017-10-10T13:19:00.000Z', '9:30am', '2013-06-19T00:00:00.000Z', '23:00-5:00', '2017-08-23T20:00:00.000Z', '8.30pm', '20:37', '20:30', '2016-09-08T12:40:00.000Z', '2017-06-27T21:04:00.000Z']
'''


