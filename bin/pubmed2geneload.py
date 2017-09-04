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
bcpon = 1        # bcp into the database?  default is yes.

refAssocKey = None
mgiTypeKey = 2	# marker
refAssocTypeKey = 1018 # General
createdByKey = 1571
loaddate = loadlib.loaddate

# {pmID:[[mgiID, jNum, refsKey], ...], ...}
dbPmToMgiDict = {}

# {egID:[pmID, ...], ...}
inputEgToPmDict = {}

# {egID:[ [m1,symbols1, markerKey] ...], ...}
dbEgToMarker = {}

# curation log
logCur = ''

# diagnostic log
logDiag = ''

# bcp file
fpBcp = ''



# for reporting to the curatiion log
inputPmIdNotInMgiList = []
inputPmIdMultiEgList = []
egIdMultiGenesList = []

def init():
    global dbPmToMgiDict, inputEgToPmDict, dbEgToMarker
    global refAssocKey, fpBcp

    # curation log
    logCur = open (os.environ['LOG_CUR'], 'a')

    # diagnostic log
    logDiag = open (os.environ['LOG_DIAG'], 'a')
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
    for r in results:
	pmID = r['pubmedid']
	egID = r['geneid']
 	if egID not in inputEgToPmDict:
	    inputEgToPmDict[egID] = []
	inputEgToPmDict[egID].append(pmID)

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
    for egId in inputEgToPmDict:
	pmList = inputEgToPmDict[egId]
	if len(pmList) > 1:
	    inputPmIdMultiEgList.append('%s%s%s%s' % (egId, TAB, string.join(pmList), CRT))
	    continue
	pmId = pmList[0]  # here, we know we only have one
	if pmId not in dbPmToMgiDict:
	    inputPmIdNotInMgiList.append('%s%s%s%s' % (egId, TAB, pmId, CRT))
	elif len(dbEgToMarker[egId]) > 1:
	    markerList = dbEgToMarker[egId]
	    reportList = []
	    for mId, symbol in markerList:
		reportList.append('%s|%s' % (mId, symbol))
	    egIdMultiGenesList.append('%s%s%s%s%s%s' % (egId, TAB, pmId, TAB, string.join(reportList), CRT))
	elif len(dbPmToMgiDict[pmId]) > 1:
	    print 'more than one if for %s in database %s' % (pmId, dbPmToMgiDict[pmId])
	else:
	    refInfo = dbPmToMgiDict[pmId][0]
	    print refInfo
	    refKey = refInfo[2]
	    markerInfo = dbEgToMarker[egId][0]
	    markerKey = markerInfo[2]
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
    fpBcp.close()
    if DEBUG or not bcpon:
	return 0

    db.commit()

    bcpCommand = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'

    bcpCmd = '%s %s %s %s %s %s "\|" "\\n" mgd' % \
        (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), strainTable, outputDir, strainTableBCP)
    diagFile.write('%s\n' % bcpCmd)
    os.system(bcpCmd)

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
    if logCur:
	logCur.close()
    if logDiag:
	logDiag.close()
    return 0

# Main

if init() != 0:
    closeFiles()
    sys.exit(1)

if createBCP() != 0:
    closeFiles()
    sys.exit(1)
#if bcpFiles != 0:
#    closeFiles()
#    sys.exit(1)

closeFiles()

print '%s' % mgi_utils.date()
sys.exit(0)
