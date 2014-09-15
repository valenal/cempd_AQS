from cempd_AQS import *

#pollGroups = { 'VOCS':{'43505':'ACROLEIN'},
#'HAPS':{'43218':'13BUTADIENE','45201':'BENZENE','43502':'FORMALDEHYDE','43503':'ACETALDEHYDE'},
#'SPEC':{'88403':'SO4'}, 'PM25_FRM':{'88101':'PM25_FRM'},'CO':{'42101':'CO'} ,
#'PM25_nonFRM':{'88502':'PM25_nonFRM'},'NOX':{'42603':'NOX'}}

pollGroups = {'SO2':{'42401':'SO2'}}

states = ['45']
counties = ['019','015']
years = [2011]
tmean = 'hourly'

for year in years:
    for polls in pollGroups.keys():
        outFile = DL_unzip(polls,year,tmean)
        nex = AQSdat(outFile,polls,year,tmean)


        for pollCode in pollGroups[polls].keys():
            pollName = pollGroups[polls][pollCode]
            nexSts = nex.keepStates(states).keepStates(counties,'county')
            spc = nexSts.getPoll(pollCode,pollName,toUgm3=True)
            spc.writeAMETRDY(createLocsFile=True)
