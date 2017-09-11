#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       Create pubmed to gene associations in the database
#
# Usage: pubmed2geneload.py
# Env Vars:
#	 1. 
#
# Inputs:
#	1. RADAR database DP_EntrezGene_PubMed table
#	2. Configuration (see pubmed2geneload.py
#
# Outputs:
#	 1. bcp file
#	 2. logs
#	 3. reports
# 
# Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
#  Assumes:  Nothing
#
#  Notes:  None
#
###########################################################################

import sys
import os
import mgi_utils
import string
import db
import loadlib

print '%s' % mgi_utils.date()

TAB = '\t'
CRT = '\n'

DEBUG = 0	 # if 0, not in debug mode
bcpon =  1       # bcp into the database?  default is yes.

refAssocKey = None
mgiTypeKey = 2	# marker
refAssocTypeKey = 1018 # General
createdByKey = 1571
loaddate = loadlib.loaddate

# {pmID:[[mgiID, jNum, refsKey], ...], ...}
dbPmToMgiDict = {}

# {egID:[pmID, ...], ...}
#inputEgToPmDict = {}
# pmID:[egID, ...], ...}
inputPmToEgDict = {}

# {egID:[ [m1,symbols1, markerKey] ...], ...}
dbEgToMarker = {}

# curation log
fpLogCur = ''

# curated reference list
#{ refsKey: markerKey, ...}
curRefDict = {}

# diagnostic log
fpLogDiag = ''

# Output directory and table name for bcp
outputDir = os.environ['OUTPUTDIR']
refTable = 'MGI_Reference_Assoc'

# bcp file
bcpFile =  refTable + '.bcp'
fpBcp = ''

# Total eg ID/reference pairs in DP_EntrezGene_PubMed
totalAssoc = 0

# Total gene/reference associations deleted from db
totalDeleted = 0

# Total gene/reference associations added to the databaase
totalAdded = 0

# total assoc already in db
totalAssocInDb = 0

# for reporting to the curation log
inputPmIdNotInMgiList = [] # 1 PM ID not in the database
inputPmIdMultiEgList = [] # 2 Reference Associated with > egID in input
inputEgIdNotInMgiList =  [] # 3 EG ID not in the database  
egIdMultiGenesList = [] # 4 EG ID associated with < 1 marker in database

def init():
    #global dbPmToMgiDict, inputEgToPmDict, dbEgToMarker
    global dbPmToMgiDict, inputPmToEgDict, dbEgToMarker, curRefDict
    global refAssocKey, fpBcp, fpLogCur, fpLogDiag, totalAssoc

    # curation log
    fpLogCur = open (os.environ['LOG_CUR'], 'a')

    # diagnostic log
    fpLogDiag = open (os.environ['LOG_DIAG'], 'a')

    # bcp file
    fpBcp = open(os.environ['BCP_FILE'], 'w')


    results = db.sql('select max(_Assoc_key) + 1 as maxKey from MGI_Reference_Assoc', 'auto')
    refAssocKey = results[0]['maxKey']

    # -- get all refs with pubmed IDs
    # -- use to determine if PM id not in db
    # -- used in report of multi egIDs assoc with pmid in file
    results = db.sql('''select a1._Object_key as refsKey, a1.accid as mgiID, 
	    a2.accid as jNum, a3.accid as pmID
	from ACC_Accession a1, ACC_Accession a2, ACC_Accession a3
	where a1._MGIType_key = 1
	and a1._LogicalDB_key = 1
	and a1.preferred = 1
	and a1.prefixPart = 'MGI:'
	and a1._Object_key = a2._Object_key 
	and a2._MGIType_key = 1
	and a2._LogicalDB_key = 1
	and a2.preferred = 1
	and a2.prefixPart = 'J:'
	and a1._Object_key = a3._Object_key
	and a3._MGIType_key = 1
	and a3._LogicalDB_key = 29
	and a3.preferred = 1''', 'auto')
    for r in results:
	pmID = r['pmID']
	mgiID = r['mgiID']
	jNum = r['jNum']
	refsKey = r['refsKey']
	if pmID not in dbPmToMgiDict:
	    dbPmToMgiDict[pmID] = []
	dbPmToMgiDict[pmID].append([mgiID, jNum, refsKey])
    # -- get mouse gene/ref pairs
    # -- load into dict, use to iterate through pubmed/id gene pairs
    # -- and report genes assoc with >1 pmid at EG
    results = db.sql('''select geneid, pubmedid
	from DP_EntrezGene_PubMed
	where taxid = 10090''', 'auto')
    totalAssoc = len(results)
    for r in results:
	pmID = r['pubmedid']
	egID = r['geneid']
 	#if egID not in inputEgToPmDict:
	#    inputEgToPmDict[egID] = []
	#inputEgToPmDict[egID].append(pmID)
	if pmID  not in inputPmToEgDict:
	    inputPmToEgDict[pmID] = []
	inputPmToEgDict[pmID].append(egID)
	

    # -- use to determine if egID assoc >1 marker
    # -- report EG ID and gene mgiID and marker symbol
    results = db.sql('''select a1.accid as egID, m._Marker_key, m.symbol, 
	    a2.accid as markerID
	from ACC_Accession a1, MRK_Marker m, ACC_Accession a2
	where a1._MGIType_key = 2
	and a1._LogicalDB_key = 55
	and a1.preferred = 1
	and a1._Object_key = m._Marker_key
	and a1._Object_key = a2._Object_key
	and a2._MGIType_key = 2
	and a2._LogicalDB_key = 1
	and a2.preferred = 1
	and a2.prefixPart = 'MGI:'
	order by a1.accid''', 'auto')
    for r in results:
	egID = r['egID']
	markerID = r['markerID']
	symbol = r['symbol']
	markerKey = r['_Marker_key']
        if egID not in dbEgToMarker:
	    dbEgToMarker[egID] = []
	dbEgToMarker[egID].append([markerID, symbol, markerKey])
    # -- used to determine if there is an existing ref/egid curated
    # -- association in the database
    # -- 1571 - exclude this loads user key as it hasn't been deleted yet.
    results = db.sql('''select distinct ra._Refs_key, a._Object_key as markerKey, a.accid as egId
	from MGI_Reference_Assoc ra, ACC_Accession a, MGI_User u
	where ra._MGIType_key = 2
	and ra._CreatedBy_key != 1571
	and ra._CreatedBy_key = u._User_key
	and ra._Object_key = a._Object_key 
	and a._MGIType_key = 2
	and a._LOgicalDB_key = 55 --entrezgene
	and a.preferred = 1''', 'auto')

    for r in results:
	refsKey = r['_Refs_key']
	egID = r['egId']
	markerKey = r['markerKey']
	if refsKey not in curRefDict:
	    curRefDict[refsKey] = []
	curRefDict[refsKey].append(markerKey)
    return 0

#
# Purpose: create reference association bcp file, produce error reports
# Returns: 0
# Assumes: input/output files exist and have been opened
# Effects: writes to the file system
# Throws: Nothing
#

def createBCP():
    global refAssocKey
    global inputPmIdNotInMgiList, inputPmIdMultiEgList, inputEgIdNotInMgiList
    global egIdMultiGenesList, totalAssocInDb, totalAdded, totalDeleted

    #for egID in inputEgToPmDict:
        #pmList = inputEgToPmDict[egID]
    for pmID in inputPmToEgDict:
	egList = inputPmToEgDict[pmID]
	#if len(pmList) > 1: # 2
        #    inputPmIdMultiEgList.append('%s%s%s%s' % (egID, TAB, string.join(pmList), CRT))
        #    continue

	if len(egList) > 1: # 2
	    inputPmIdMultiEgList.append('%s%s%s' % (pmID, TAB, string.join(egList)))
	    continue
	egID = egList[0]  # here, we know we only have one
	if pmID not in dbPmToMgiDict: # 1
	    inputPmIdNotInMgiList.append('%s%s%s' % (egID, TAB, pmID))
	    continue
	if egID not in dbEgToMarker: # 3
	    inputEgIdNotInMgiList.append('%s%s%s' % (egID, TAB, pmID))
            continue
	if len(dbEgToMarker[egID]) > 1: # 4
	    markerList = dbEgToMarker[egID]
	    reportList = []
	    for l in markerList:
		mID = l[0]
		symbol = l[1]
		reportList.append('%s|%s' % (mID, symbol))
	    egIdMultiGenesList.append('%s%s%s%s%s' % (egID, TAB, pmID, TAB, string.join(reportList)))
	    continue
	if len(dbPmToMgiDict[pmID]) > 1:
	    # Just curious so checking  - this will go to diag log
	    print 'more than one reference object for %s in database %s' % (pmID, dbPmToMgiDict[pmID])
	    continue

	refInfo = dbPmToMgiDict[pmID][0]
	#print refInfo
	refKey = refInfo[2]
	markerInfo = dbEgToMarker[egID][0] # [mgiID, jNum, refKey]
	markerKey = markerInfo[2]
	# skip if the reference/marker assoc already curated in the db
	if refKey in curRefDict and markerKey in curRefDict[refKey]:
	    totalAssocInDb += 1
	    continue
	totalAdded += 1
	# START HERE 2 call updateGoStatus with refKey
	fpBcp.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    % (refAssocKey, refKey, markerKey, mgiTypeKey, refAssocTypeKey, \
	    createdByKey, createdByKey, loaddate, loaddate))
	refAssocKey = refAssocKey + 1
    return 0

# Purpose:  BCPs the data into the database
# Returns:  nothing
# Assumes:  nothing
# Effects:  BCPs the data into the database
# Throws:   nothing

def bcpFiles():
    global totalDeleted
    fpBcp.close()
    if DEBUG or not bcpon:
	return 0

    db.commit()

    # delete from MGI_Reference_Assoc
    db.sql('''select _Assoc_key
	into temporary table toDelete
	from MGI_Reference_Assoc
	where _CreatedBy_key = 1571''', None)

    db.sql('''create index idx1 on toDelete(_Assoc_key)''', None)

    results = db.sql('''select count(*) as deleteCt
	from toDelete''', 'auto')

    for r in results:
	totalDeleted = r['deleteCt']
    print 'totalDeleted in bcpFiles(): %s' % totalDeleted
    db.sql('''delete from MGI_Reference_Assoc a
	using toDelete d
	where a._Assoc_key = d._Assoc_key''', None)

    db.commit()
    # add new associations
    bcpCommand = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'

    bcpCmd = '%s %s %s %s %s %s "\|" "\\n" mgd' % \
        (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), refTable, outputDir, bcpFile)
    print bcpCmd
    fpLogDiag.write('%s\n' % bcpCmd)
    print 'executing bcp'
    os.system(bcpCmd)

    return 0

# START HERE 1
#def updateGoStatus(refKey):
    # for all reference associations we are creating:
    # create lookup to check existing GO group status - 
    # 	if it is 'rejected' or 'discard', report
    # if current status is not 'Fully Coded' (or rejected or discard 
    # which is reported above), set status to 'Indexed' and generate Jnum
    # call the API to update the GO group status and  create Jnum 
#
# Purpose: write all errors to curation log
# Returns: Nothing
# Assumes: Nothing
# Effects: writes to curation log
# Throws: Nothing
#
def writeCuratorLog():

    #inputPmIdNotInMgiList = [] # 1 PM ID not in the database 
	# egID, TAB, pmID, CRT
    #inputPmIdMultiEgList = [] # 2 Reference Associated with > egID in input 
	# pmID, TAB, string.join(egList), CRT
    #inputEgIdNotInMgiList =  [] # 3 EG ID not in the database 
	# egID, TAB, pmID, CRT
    #egIdMultiGenesList = [] # 4 EG ID associated with < 1 marker in database 
	# egID, TAB, pmID, TAB, string.join(reportList)
    print 'totalDeleted in writeCuratorLog: %s' % totalDeleted
    fpLogCur.write(CRT + CRT + 'Total PubMed/EG ID associations from EntrezGene: %s' % totalAssoc)
    fpLogCur.write(CRT + CRT + 'Total PubMed/EG ID associations already in Database: %s' % totalAssocInDb)
    fpLogCur.write(CRT + CRT + 'Total PubMed/EG ID associations deleted from the Database: %s' % totalDeleted)
    fpLogCur.write(CRT + CRT + 'Total PubMed/EG ID added to the Database: %s' % totalAdded)
    if len(inputPmIdNotInMgiList):
	fpLogCur.write(CRT + CRT + string.center(
	    'PM IDs not in the Database',60) + CRT)
	fpLogCur.write('%-12s  %-20s%s' %
             ('EG ID','PM ID', CRT))
	fpLogCur.write(string.join(inputPmIdNotInMgiList, CRT))
	fpLogCur.write(CRT + 'Total: %s' % len(inputPmIdNotInMgiList))
    if len(inputPmIdMultiEgList):
	fpLogCur.write(CRT + CRT + string.center(
            'PM IDs  Associated with > 1 egID in Input',60) + CRT)
	fpLogCur.write('%-12s  %-20s%s' %
            ('PM ID','EG IDs', CRT))
	fpLogCur.write(string.join(inputPmIdMultiEgList, CRT))
	fpLogCur.write(CRT + 'Total: %s' % len(inputPmIdMultiEgList))
    if len(inputEgIdNotInMgiList):
	fpLogCur.write(CRT + CRT + string.center(
            'EG IDs not in the Database',60) + CRT)
	fpLogCur.write('%-12s  %-20s%s' %
             ('EG ID','PM ID', CRT))
	fpLogCur.write(string.join(inputEgIdNotInMgiList, CRT))
	fpLogCur.write(CRT + 'Total: %s' % len(inputEgIdNotInMgiList))
    if len(egIdMultiGenesList):
	fpLogCur.write(CRT + CRT + string.center(
            'EG IDs associated with < 1 marker in the Database',60) + CRT)
	fpLogCur.write('%-12s  %-20s  %-20s%s' %
             ('EG ID','PM ID', 'Markers', CRT))
	fpLogCur.write(string.join(egIdMultiGenesList, CRT))
	fpLogCur.write(CRT + 'Total: %s' % len(egIdMultiGenesList))
    return 0
#
# Purpose: Close files.
# Returns: 0
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def closeFiles():
    if fpBcp:
	fpBcp.close()
    if fpLogCur:
	fpLogCur.close()
    if fpLogDiag:
	fpLogDiag.close()
    return 0

# Main

if init() != 0:
    print 'Initialization failed'
    closeFiles()
    sys.exit(1)

if createBCP() != 0:
    print 'Creating BCP Files failed'
    closeFiles()
    sys.exit(1)

if bcpFiles() != 0:
    print 'BCP failed'
    closeFiles()
    sys.exit(1)

if writeCuratorLog() != 0:
    print 'Writing to curator log failed'
    closeFiles()
    sys.exit(1)

closeFiles()

print '%s' % mgi_utils.date()
sys.exit(0)
