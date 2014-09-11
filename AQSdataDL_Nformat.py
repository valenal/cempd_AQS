import os,time,pdb
import urllib
import urllib2
import zipfile
import pandas as pd
import numpy as np
import datetime as dt
from IPython import embed

#######################################################################
# Conversion factors are ppb -> ug/m3
# Conversion factors are calculated by the equation
# Density of Air at 25C (1.18535 Standard) * (mw conc / mw dry Air (28.9652))
# http://www.engineeringtoolbox.com/air-density-specific-weight-d_600.html
#######################################################################
AQSpollsDict ={'42101':dict(name="CO", convert=1.14626, nCarbon=0.0),
'42602':dict(name="NO2", convert=1.88269, nCarbon=0.0),
'42603':dict(name="NOX", convert=1.88288, nCarbon=0.0),
'42600':dict(name="NOy", convert=1.88288, nCarbon=0.0),
'42401':dict(name="SO2", convert=2.860, nCarbon=0.0),
'88101':dict(name="PM25", convert=-999, nCarbon=0.0),
'81102':dict(name="PM10", convert=-999, nCarbon=0.0),
'45201':dict(name="BENZENE", convert=3.1965, nCarbon=6.0),
'88305':dict(name="OC25", convert=-999, nCarbon=0.0),
'88307':dict(name="EC25", convert=-999, nCarbon=0.0),
'88403':dict(name="SO4", convert=-999, nCarbon=0.0),
'88501':dict(name="PM25", convert=-999, nCarbon=0.0),
'88502':dict(name="PM25", convert=-999, nCarbon=0.0),
'44201':dict(name="O3", convert=1.9643158, nCarbon=0.0),
'43502':dict(name="FORMALDEHYDE", convert=1.22892, nCarbon=1.0),
'43218':dict(name="1,3-BUTADIENE", convert=2.213604, nCarbon=4.0),
'43505':dict(name="ACROLEIN", convert=2.294157, nCarbon=3.0),
'43509':dict(name="ACROLEIN", convert=2.294157, nCarbon=3.0),
'43503':dict(name="ACETALDEHYDE", convert=1.80277, convert0=1.96752, nCarbon=2.0)}
       
units = {"Parts per million":"ppm",
"Parts per billion Carbon":"ppbC",
"Micrograms/cubic meter (LC)":"ug/m3",
"Micrograms/cubic meter (25 C)":"ug/m3",
"Nanograms/cubic meter (25 C)":'ng/m3'}

basePath = os.path.dirname(__file__)


### FIX THIS
def getAdjStatesDict(getCode=True):
    f1 = open('%s/adjacent_states.txt' % basePath)
    f2 = open('%s/statesCodeAbbName.txt' % basePath)
    
    #get State Code
    if getCode:
        iCode = 1
    #get State Name
    else:
        iCode = 2 
        
    #create dict that map abbreviations to State code or name
    abbDict = dict([ (i.rstrip('\n').split(',')[0],i.rstrip('\n').split(',')[iCode])  for i in f2.readlines()])    
    #create dict with adjacent states
    outDict = dict([ (i.rstrip('\n').split(',')[0],i.rstrip('\n').split(',')[1:])  for i in f1.readlines()])

    f1.close()
    f2.close()

    for key,vals in outDict.items():
        outDict[key] = [ abbDict[i] for i in vals]
              
    return outDict


def unzip(inPath,tave,spc,code,yr,outPath):
    zfile = zipfile.ZipFile(inPath)
    for name in [i for i in zfile.namelist() if "csv" in i]:
        (dirname, filename) = os.path.split(name)
        print "Decompressing " + filename 
        if filename == '':
            # directory
            if not os.path.exists(dirname):
                os.mkdir(dirname)
        else:
            # file
            fd = open(name, 'w')
            fd.write(zfile.read(name))
            fd.close()
            os.rename(name,outPath)
    zfile.close()
    os.remove(inPath)

def DL_unzip(spc,yr,tave):

    spcCode = {'CO':['42101'],
    'NOX':['42603'], 
    'NO2':['42602'], 
    'O3':['44201'], 
    'PM25_FRM':['88101'],
    'PM25_NONFRM':['88502'],
    'PM10':['81102'], 
    'SO2':['42401'],
    'VOCS':['VOCS'],
    'SPEC':['SPEC'],
    'HAPS':['HAPS'],
    'LEAD':['LEAD']}

    if spc.upper() == 'ALL':
        spcs = spcCode.keys()
    else:
        spcs= [spc.upper()]

    for each_spcs in spcs:
        for each_spc in spcCode[each_spcs]:

            dlFile = 'http://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/%s_%s_%s.zip' % (tave,each_spc,yr)
            zipFile = '%s/AQSFiles/%s_%s_%s.zip' % (os.getcwd(),tave,each_spc,yr)
            if not os.path.exists(os.path.dirname(zipFile)):
                os.makedirs(os.path.dirname(zipFile))
            
            if each_spcs == each_spc:
                unzipFile = '%s/AQSFiles/%s_%s_%s.csv' % (os.getcwd(),tave,each_spc,yr)
            else:
                unzipFile = '%s/AQSFiles/%s_%s_%s_%s.csv' % (os.getcwd(),tave,each_spcs,each_spc,yr)
            
            if not os.path.isfile(unzipFile):
                print 'Downloading '+ each_spc +' at '+ dlFile
                resp = urllib2.urlopen(dlFile)
                urllib.urlretrieve(dlFile,zipFile)
                urllib.urlcleanup()
                print 'Unzipping '+ each_spc +' at '+ zipFile
                unzip(zipFile,tave,each_spcs,each_spc,yr,unzipFile)
            else:
                print 'File %s already downloaded on %s' % (unzipFile,time.ctime(os.path.getmtime(unzipFile)))
    return unzipFile


class AQSdat:

    def __init__(self,inPath,spc,yr,tave):
        self.pollutant = spc 
        self.year = yr 
        self.timeAve = tave 
        self.df = pd.read_csv(inPath)
        self.df['site_id'] = ['%02i%03i%04i' % (i[0],i[1],i[2]) for i in self.df.loc[:,['State Code','County Code','Site Num']].values]
            
    def keepStates(self,states,stateType='name'):
        #need to make it check for type of states inserted
        if stateType == 'Name':
            cond  = self.df['State Name'].str.contains('|'.join([i.capitalize() for i in states]))
        else:
            cond  = np.in1d(self.df['State Code'],np.array(states))
        self.df = self.df.drop(self.df[~cond].index.values)

    def getUnqPolls(self):
        return self.df.loc[:,['Parameter Name','Parameter Code']].drop_duplicates().values 

    def getPoll(self,paramCode,toUgm3=False,toPpm=False,toPpb=False):
        outDF = self.df[self.df['Parameter Code'] == int(paramCode)]
        
        try:
            fPoll = AQSpollsDict[paramCode]
        except KeyError as (strerror):
            raise KeyError("Cannot process AQS param code {0}".format(strerror))

        if self.timeAve == 'daily':
            datCol = 'Arithmetic Mean'
        elif self.timeAve == 'hourly':
            datCol = 'Sample Measurement'
        
        ###################################################
        # Create array of unit names based on lookup table
        ###################################################
        oldUnits = outDF['Units of Measure'].values.astype(str)
        newUnits = np.copy(oldUnits)
        for k, v in units.iteritems(): newUnits[oldUnits==k] = v
        outDF['Units of Measure'] = newUnits
        #Old way didnt work
        #outDF = outDF.replace({'Units of Measure':units})
        fUnitsAll = outDF['Units of Measure'].values
        unitsDiff = np.setdiff1d(fUnitsAll,np.array(units.values()) )
        if unitsDiff:
            raise Exception("No units found for %s" % unitsDiff)
        concAll = outDF[datCol].values
 
        ############################
        # Convert units as required
        # ppbC = ppb * # of Carbon atoms
        # ppm = ppb / 1000.
        ############################
        if toUgm3: # convert to ug/m3
            # convert ppm units
            b = fUnitsAll == "ppm"
            concAll[b] = concAll[b] * 1000. * fPoll['convert']

            # convert pptm units
            b = fUnitsAll == "pptm"
            concAll[b] = concAll[b] * 100. * fPoll['convert']

            # convert pphm units
            b = fUnitsAll == "pphm"
            concAll[b] = concAll[b] * 10. * fPoll['convert']

            # convert ppb units
            b = fUnitsAll == "ppb"
            concAll[b] = concAll[b] * fPoll['convert']

            # convert ppbC units
            b = fUnitsAll == "ppbC"
            concAll[b] = (concAll[b] / fPoll['nCarbon']) * fPoll['convert']

            # convert mg/m3 units
            b = fUnitsAll == "mg/m3"
            concAll[b] = concAll[b] * 1000.

            # Set units flag
            fUnits = "ug/m3"
        
        if toPpm: # convert to ppm
            # convert ppb units
            b = fUnitsAll == "ppb"
            concAll[b] = concAll[b] / 1000.

            # convert pptm units
            b = fUnitsAll == "pptm"
            concAll[b] = concAll[b] / 100.

            # convert pphm units
            b = fUnitsAll == "pphm"
            concAll[b] = concAll[b] / 10.

            # convert ppbC units
            b = fUnitsAll == "ppbC"
            concAll[b] = (concAll[b] / fPoll['nCarbon']) / 1000.

            # convert ug/m3 units
            b = fUnitsAll == "ug/m3"
            concAll[b] = concAll[b] / fPoll['convert'] / 1000.

            # convert mg/m3 units
            b = fUnitsAll == "mg/m3"
            concAll[b] = concAll[b] / fPoll['convert']

            fUnits = "ppm"

        if toPpb: # convert to ppb
            # convert ppm units
            b = fUnitsAll == "ppm"
            concAll[b] = concAll[b] * 1000.


            # convert pptm units
            b = fUnitsAll == "pptm"
            concAll[b] = concAll[b] * 100.

            # convert pphm units
            b = fUnitsAll == "pphm"
            concAll[b] = concAll[b] * 10.

            # convert ppbC units
            b = fUnitsAll == "ppbC"
            concAll[b] = concAll[b] / fPoll['nCarbon']

            # convert mg/m3 units
            b = fUnitsAll == "mg/m3"
            concAll[b] = 1000 * concAll[b] / fPoll['convert']

            # convert ug/m3 units
            b = fUnitsAll == "ug/m3"
            concAll[b] = concAll[b] / fPoll['convert']

            fUnits = "ppb"
        
        outDF['Units of Measure'] = fUnits
        #outDF[datCol] = concAll
        
        if outDF.shape[0] == 0:
            print "Warning: Empty DataFrame for %s" % fPoll['name']
        return outDF

def getUnqSites(inDF):
    return inDF['site_id'].drop_duplicates().values 

def writeEXT(inDF,spc,tmean,outF):
    if tmean == 'daily':
        inDF['dateon'] = pd.to_datetime(inDF['Date Local'])
        inDF=inDF.rename(columns = {'Arithmetic Mean':spc})        
    elif tmean == 'hourly':
        inDF['dateon'] = pd.to_datetime(inDF['Date Local']+' '+ inDF['Time Local'])
        inDF=inDF.rename(columns = {'Sample Measurement':spc})        
    inDF['dateoff'] =(inDF['dateon'] + dt.timedelta(hours=23,minutes=59))
    print "Creating %s File" % outF
    inDF.loc[:,['site_id','dateon','dateoff',spc]].to_csv(outF,index=False)

def writeSTOK(inDF,spc,tmean,outF):
    if tmean == 'daily':
        inDF['dt.date'] = pd.to_datetime(inDF['Date Local'])
        inDF=inDF.rename(columns = {'Arithmetic Mean':spc})  
    elif tmean == 'hourly':
        inDF['dt.date'] = pd.to_datetime(inDF['Date Local']+' '+ inDF['Time Local'])
        inDF=inDF.rename(columns = {'Sample Measurement':spc})          

    inDF['year']   = inDF['dt.date'].map(lambda x: x.year)
    inDF['month']   = inDF['dt.date'].map(lambda x: x.month)
    inDF['day']   = inDF['dt.date'].map(lambda x: x.day)
    inDF['hour']   = inDF['dt.date'].map(lambda x: x.hour)

    # remove POC duplicates
    inDF = inDF.groupby(['site_id','Longitude','Latitude','year','month','day','hour'],as_index=False).mean()
    #remove SAME ID slighlty different coordinate problem
    unq_sites = inDF.loc[:,['site_id','Longitude','Latitude']].drop_duplicates()
    dupli_sites = unq_sites[unq_sites.duplicated('site_id')]

    if dupli_sites.shape[0] != 0:
        for dupli_site in dupli_sites.values():
            inDF[inDF['site_id']== dupli_site_id[0]]['Longitude'] == dupli_site_id[1] 
            inDF[inDF['site_id']== dupli_site_id[0]]['Latitude']  == dupli_site_id[2] 

    print "Creating %s File" % outF
    inDF.loc[:,['site_id','Longitude','Latitude','year','month','day','hour',spc]].to_csv(outF,index=False,header=False)

if __name__ == "__main__":

######################## EXAMPLE RUN ##############################
    if (FALSE):
        pollGroups = { 'VOCS':{'43505':'ACROLEIN'},
        'HAPS':{'43218':'13BUTADIENE','45201':'BENZENE','43502':'FORMALDEHYDE','43503':'ACETALDEHYDE'},
        'SPEC':{'88403':'SO4'}, 'PM25_FRM':{'88101':'PM25_FRM'},'CO':{'42101':'CO'} ,
        'PM25_nonFRM':{'88502':'PM25_nonFRM'},'NOX':{'42603':'NOX'}}

        states = [37,51,21,47,13,45]
        years = [2011,2012,2013,2014]
        tmean = 'hourly'

        for year in years:
            for polls in pollGroups.keys():
                outFile = DL_unzip(polls,year,tmean)
                nex = AQSdat(outFile,polls,year,tmean)


                for pollCode in pollGroups[polls].keys():
                    pollName = pollGroups[polls][pollCode]
                    nex.keepStates(states,'Code')
                    spc = nex.getPoll(pollCode,toUgm3=True)
                    writeSTOK(spc,pollName,tmean,'./%s_%s.csv' % (pollName,year))
####################################################################################
