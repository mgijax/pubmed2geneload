
##########################################################################
#
# Purpose:
#       Create pubmed to gene associations in the database
#
# Usage: pubmed2geneload.py
#
# Env Vars:
#	 1. OUTPUTDIR - where to write the bcp file
#	 2. LOG_CUR - curation log for reporting
#	 3. LOG_DIAG - diagnostic log 
#	 4. BCP_FILE -  full path to the bcp file
#	 5. PG_DBUTILS - for the bcp utility script
#	 6. JAVA_API_URL - access to the java api for status updates
#	 7. JAVA_API_TOKEN - security token for java api updates
#	 8. UPDATE_BATCH - size of the status update batches
#
# Inputs:
#	1. RADAR database DP_EntrezGene_PubMed table
#	2. Configuration (see pubmed2geneload.config)
#
# Outputs:
#	 1. bcp file
#	 2. logs
# 
# Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
#  Assumes:  Nothing
#
#  Notes: This script uses the Java API (mgd_java_api) to do status updates
#
# History
#	2/9/2018	sc
#	- TR12760 - create marker/reference associations and update GO status
#		only if reference has JNumber
###########################################################################

import sys
import os
import mgi_utils
import string
import db
import loadlib
import runCommand
        
TAB = '\t'
CRT = '\n'

DEBUG = 0	 # if 0, not in debug mode
bcpon = 1        # bcp into the database?  default is yes.

# for creating API URL to update status
PUT = 'curl -X PUT "'
JAVA_API_URL = os.environ['JAVA_API_URL']
JAVA_API_TOKEN = os.environ['JAVA_API_TOKEN']
STATUS_URL = 'reference/statusUpdate?accid=%s&group=GO&status=Indexed" -H "accept: application/json" -H  "api_access_token: ' + JAVA_API_TOKEN + '" -H  "username: pm2geneload"'

# the full URL - just plug in comma delim (no space) list of reference MGI ID
FULL_API_URL = PUT + JAVA_API_URL + STATUS_URL

# max number of status updates in a batch
UPDATE_BATCH = int(os.environ['UPDATE_BATCH'])

# next available MGI_Reference_Assoc primary key
refAssocKey = None
mgiTypeKey = 2	# marker
refAssocTypeKey = 1018 # General
createdByKey = 1571 # pubmed2geneload
loaddate = loadlib.loaddate # today's date

# {pmID:[[mgiID, refsKey], ...], ...}
dbPmToMgiDict = {}

# pmID:[egID, ...], ...}
inputPmToEgDict = {}

# {egID:[ [m1,symbols1, markerKey] ...], ...}
dbEgToMarkerDict = {}

# {refID:[isDiscard, _Status_key], ...}
dbRefIdToStatusDict = {}

# curation log
fpLogCur = ''

# curated reference list
#{ refsKey: markerKey, ...}
curRefDict = {}

# _refs_key,_object_key must be unique
# use this list to check for duplicates
# duplicates will be skipped
refMarkerList = []

# list of status updates
statusUpdateList = []

# diagnostic log
fpLogDiag = ''

# Output directory and table name for bcp
outputDir = os.environ['OUTPUTDIR']
refTable = 'MGI_Reference_Assoc'

# bcp file
bcpFile =  refTable + '.bcp'
fpBcp = ''

# Total eg ID/reference pairs in DP_EntrezGene_PubMed
# Total gene/reference associations deleted from db
totalDeleted = 0

# Total gene/reference associations added to the databaase
totalAdded = 0

# total mouse associations in input
totalAssocInput = 0

# total assoc already in db
totalAssocInDb = 0

# references that will be added to the database and potentially need
# status updates
refList = []

# list of references that actually need status updates
updateStatusList = []
numUpdates = 0

# lists of discrepancies for reporting to the curation log
inputPmIdNotInMgiList = [] # 1 PM ID not in the database
inputPmIdMultiEgList = [] # 2 Reference Associated with > 15 egID in input
inputEgIdNotInMgiList =  [] # 3 EG ID not in the database  
egIdMultiGenesList = [] # 4 EG ID associated with > 1 marker in database
dOrRStatusList = [] # 5 current status is discard or rejected

#
# Purpose: initialize lookups, next primary key, file descriptors
# Returns: 0
# Assumes: input/output files exist and have been opened
# Effects: writes to the file system
# Throws: Nothing
#
def init():
    global dbPmToMgiDict, inputPmToEgDict, dbEgToMarkerDict, curRefDict
    global refAssocKey, fpBcp, fpLogCur, fpLogDiag, totalAssocInput
    global dbRefIdToStatusDict

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # curation log
    fpLogCur = open (os.environ['LOG_CUR'], 'a')

    # diagnostic log
    fpLogDiag = open (os.environ['LOG_DIAG'], 'a')

    # bcp file
    fpBcp = open(os.environ['BCP_FILE'], 'w')


    results = db.sql(''' select nextval('mgi_reference_assoc_seq') as maxKey ''', 'auto')
    refAssocKey = results[0]['maxKey']

    # -- get all refs with pubmed IDs AND JNumbers
    # -- use to determine if pm ID associated w/reference that has a JNum
    # -- used in report of multi egIDs assoc with pmid in file
    results = db.sql('''select a1._Object_key as refsKey, a1.accid as mgiID,
             a3.accid as pmID
        from ACC_Accession a1, ACC_Accession a2, ACC_Accession a3
        where a1._MGIType_key = 1
        and a1._LogicalDB_key = 1
        and a1.preferred = 1
        and a1.prefixPart = 'MGI:'
        and a1._Object_key = a2._Object_key
        and a2._MGIType_key = 1
        and a2._LogicalDB_key = 1
        and a2.prefixPart = 'J:'
        and a2.preferred = 1
        and a1._Object_key = a3._Object_key
        and a3._MGIType_key = 1
        and a3._LogicalDB_key = 29
        and a3.preferred = 1''', 'auto')
    for r in results:
        pmID = r['pmID']
        mgiID = r['mgiID']
        refsKey = r['refsKey']
        if pmID not in dbPmToMgiDict:
            dbPmToMgiDict[pmID] = []
        dbPmToMgiDict[pmID].append([mgiID, refsKey])
    # -- get mouse gene/ref pairs
    # -- load into dict, use to iterate through pubmed/id gene pairs
    # -- and report genes assoc with >1 pmid at EG
    results = db.sql('''select geneid, pubmedid
        from DP_EntrezGene_PubMed
        where taxid = 10090''', 'auto')
    totalAssocInput = len(results)
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
        if egID not in dbEgToMarkerDict:
            dbEgToMarkerDict[egID] = []
        listToAppend = [markerID, symbol, markerKey]
        if listToAppend not in dbEgToMarkerDict[egID]:
            dbEgToMarkerDict[egID].append(listToAppend)
    # -- used to determine if there is an existing ref/egid curated
    # -- association in the database
    # -- 1571 - exclude this loads user key as it hasn't been deleted yet.
    results = db.sql('''select distinct ra._Refs_key, a._Object_key as markerKey, a.accid as egId
        from MGI_Reference_Assoc ra, ACC_Accession a, MGI_User u
        where ra._MGIType_key = 2
        and ra._refassoctype_key = 1018
        and ra._CreatedBy_key != 1571
        and ra._CreatedBy_key = u._User_key
        and ra._Object_key = a._Object_key 
        and a._MGIType_key = 2
        and a._LogicalDB_key = 55 --entrezgene
        and a.preferred = 1''', 'auto')

    for r in results:
        refsKey = r['_Refs_key']
        egID = r['egId']
        markerKey = r['markerKey']
        if refsKey not in curRefDict:
            curRefDict[refsKey] = []
        curRefDict[refsKey].append(markerKey)

    # -- get status of all refs in the database
    # -- used to update GO status when refs associated with genes
    db.sql('''select _Refs_key, _Status_key
        into temporary table goCurrent
        from BIB_Workflow_Status
        where _Group_key = 31576666
        and isCurrent = 1''', None)
    db.sql('create index idx2 on goCurrent(_Refs_key)', None)
    results = db.sql('''select a.mgiid as refID, b._Refs_key, b.isDiscard, 
            s._Status_key
        from BIB_Refs b, goCurrent s, BIB_Citation_Cache a
        where b._Refs_key = s._Refs_key
        and b._Refs_key = a._Refs_key''', 'auto')
    for r in results:
        refID = r['refID']
        isDiscard = r['isDiscard']
        statusKey = r['_Status_key']
        dbRefIdToStatusDict[refID] = [isDiscard, statusKey]
    return 0

#
# Purpose: create reference association bcp file, writes discrepancies
#	to data structures by type
# Returns: 0
# Assumes: input/output files exist and have been opened
# Effects: writes to the file system
# Throws: Nothing
#
def createBCP():
    global refAssocKey, inputPmIdNotInMgiList, inputPmIdMultiEgList
    global inputEgIdNotInMgiList,egIdMultiGenesList, totalAssocInDb
    global totalAdded, totalDeleted, refList, refMarkerList
        
    for pmID in inputPmToEgDict:
        egList = inputPmToEgDict[pmID]

        if len(egList) > 15: # 2
        inputPmIdMultiEgList.append('%s%s%s' % (pmID, TAB, ''.join(egList)))
            continue
        #egID = egList[0]  # we know we have < 16
        if pmID not in dbPmToMgiDict: # 1
            inputPmIdNotInMgiList.append(pmID)
            continue
        if len(dbPmToMgiDict[pmID]) > 1:
            # Just curious so checking  - this will go to diag log
            print('more than one reference object for %s in database %s' % (pmID, dbPmToMgiDict[pmID]))
            continue

        # get the reference info first
        refInfo = dbPmToMgiDict[pmID][0] # Here we now we have one as multi's
                                         # reported above
                                         # [mgiID, refsKey]
        refID  = refInfo[0]
        refKey = refInfo[1]

        # now iterate over the genes for this reference
        for egID in egList:
            # egID not in MGI?, report
            if egID not in dbEgToMarkerDict: # 3
                inputEgIdNotInMgiList.append('%s%s%s' % (egID, TAB, pmID))
                continue
            # egID associated with > marker in MGI?, report
            if len(dbEgToMarkerDict[egID]) > 1: # 4
                markerList = dbEgToMarkerDict[egID]
                reportList = []
                for l in markerList:
                    mID = l[0]
                    symbol = l[1]
                    reportList.append('%s|%s' % (mID, symbol))
                egIdMultiGenesList.append('%s%s%s%s%s' % \
                    (egID, TAB, pmID, TAB, str.join(reportList)))
                continue

            # get the marker  key
            markerInfo = dbEgToMarkerDict[egID][0] # Here we know we have one
                                                   # as multi reported above
            
            markerKey = markerInfo[2]
            # skip if the reference/marker assoc already curated in the db
            if refKey in curRefDict and markerKey in curRefDict[refKey]:
                totalAssocInDb += 1
                continue

            # skip if reference/marker is a duplicate
            if (refKey, markerKey) in refMarkerList:
                continue
            else:
                refMarkerList.append((refKey, markerKey))

            if refID not in refList:
                refList.append(refID)
            totalAdded += 1
            fpBcp.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' % \
                   (refAssocKey, refKey, markerKey, mgiTypeKey, refAssocTypeKey,  createdByKey, createdByKey, loaddate, loaddate))
            refAssocKey = refAssocKey + 1
    return 0

#
# Purpose:  BCPs the data into the database
# Returns:  0
# Assumes:  nothing
# Effects:  BCPs the data into the database
# Throws:   nothing
#
def bcpFiles():
    global totalDeleted, refsDeletedList
    fpBcp.close()
    if DEBUG or not bcpon:
        return 0

    db.commit()

    # delete from MGI_Reference_Assoc
    print('deleting from MGI_Reference_Assoc')
    db.sql('''select _Assoc_key
        into temporary table toDelete
        from MGI_Reference_Assoc
        where _CreatedBy_key = 1571''', None)

   
    db.sql('''create index idx1 on toDelete(_Assoc_key)''', None)

    results = db.sql('''select count(*) as deleteCt
        from toDelete''', 'auto')

    for r in results:
        totalDeleted = r['deleteCt']

    db.sql('''delete from MGI_Reference_Assoc a
        using toDelete d
        where a._Assoc_key = d._Assoc_key''', None)

    db.commit()
    print('totalDeleted in bcpFiles(): %s' % totalDeleted)

    # add new associations
    bcpCommand = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'

    bcpCmd = '%s %s %s %s %s %s "\|" "\\n" mgd' % \
        (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), refTable, outputDir, bcpFile)
    print(bcpCmd)
    fpLogDiag.write('%s\n' % bcpCmd)
    print('executing bcp')
    os.system(bcpCmd)

    # update mgi_reference_assoc auto-sequence
    db.sql(''' select setval('mgi_reference_assoc_seq', (select max(_assoc_key) from MGI_Reference_Assoc)) ''', None)
    db.commit()

    return 0

#
# Purpose:  Determines references that need GO status updates; updates them
#	    in the database using Java API calls
# Returns:  0 if successful, 1 if update failed
# Assumes:  refList has been loaded
# Effects:  Updates GO status in the database
# Throws:   nothing
#
def updateGoStatus():
    global refList, dOrRStatusList, numUpdates, updateStatusList
    # for all reference associations we are creating:
    #   check existing GO group status - 
    # 	  if it is 'rejected' or is isDiscard=true for GO, report
    #     if current status not 'Fully Coded' for GO set status to 'Indexed' 
    #       and generate Jnum, if necessary, using Java API

    # get list of refIDs that actually need updating, reporting as we go
    print('updateGoStatus size of refList: %s' % len(refList))
    for refID in refList:
        if refID in dbRefIdToStatusDict:
            print('updateGoStatus: %s' % refID)
            infoList = dbRefIdToStatusDict[refID]
            isDiscard = infoList[0]
            statusKey = infoList[1]
            print('refID: %s, isDiscard: %s statusKey: %s' % (refID, isDiscard, statusKey))
            if isDiscard == 1 or statusKey == 31576672: # rejected
                # report
                d = 'False'
                s = 'Not Rejected'
                if isDiscard == 1:
                    d = 'True'
                if statusKey == 31576672:
                    s = 'Rejected'
                dOrRStatusList.append('%s%s%s%s%s' % (refID, TAB, d, TAB, s))
                print('writing %s to discardOrRejected report' % refID)
                continue
            if statusKey not in (31576673, 31576674): # Indexed, Full-coded
                print('adding %s to list to be updated' % refID)
                updateStatusList.append(refID)
        else: # this should never happen
            print('RefID not in dbRefIdToStatusDict: %s' % refID)
    # now do the updates
    numUpdates = len(updateStatusList)
    print('updateGOStatus size of updateStatusList: %s' % len(updateStatusList))
    while updateStatusList != []:
        batchToRun = ','.join(updateStatusList[0:UPDATE_BATCH])
        #print batchToRun
        del updateStatusList[0:UPDATE_BATCH]
        print('%s' % mgi_utils.date())
        print('running runCommand')
        print(FULL_API_URL % batchToRun)
        stdout, stderr, returnCode = runCommand.runCommand(FULL_API_URL % batchToRun)
        print('after runCommand stdout: %s stderr: %s returnCode: %s' % (stdout, stderr, returnCode))
        if returnCode != 0:
            return 1
    return 0
#
# Purpose: write all discrepancies to curation log
# Returns: 0
# Assumes: Nothing
# Effects: writes to curation log
# Throws: Nothing
#
def writeCuratorLog():
    fpLogCur.write(CRT + CRT + 'Total PubMed/EG ID associations from EntrezGene: %s' % totalAssocInput)
    fpLogCur.write(CRT + CRT + 'Total PubMed/EG ID associations already in Database: %s' % totalAssocInDb)
    fpLogCur.write(CRT + CRT + 'Total PubMed/EG ID associations deleted from the Database: %s' % totalDeleted)
    fpLogCur.write(CRT + CRT + 'Total PubMed/EG ID added to the Database: %s' % totalAdded)
    fpLogCur.write(CRT + CRT + 'Total References  updated to "Indexed" for GO: %s' % numUpdates)
    if len(inputPmIdNotInMgiList):
        fpLogCur.write(CRT + CRT + str.center(
            'PM IDs not in the Database or reference has no J: number',60) + CRT)
        fpLogCur.write(CRT.join(inputPmIdNotInMgiList))
        fpLogCur.write(CRT + 'Total: %s' % len(inputPmIdNotInMgiList))
    if len(inputPmIdMultiEgList):
        fpLogCur.write(CRT + CRT + str.center(
            'PM IDs  Associated with > 15 egID in Input',60) + CRT)
        fpLogCur.write('%-12s  %-20s%s' %
            ('PM ID','EG IDs', CRT))
        fpLogCur.write(CRT.join(inputPmIdMultiEgList))
        fpLogCur.write(CRT + 'Total: %s' % len(inputPmIdMultiEgList))
    if len(inputEgIdNotInMgiList):
        fpLogCur.write(CRT + CRT + str.center(
            'EG IDs not in the Database',60) + CRT)
        fpLogCur.write('%-12s  %-20s%s' %
             ('EG ID','PM ID', CRT))
        fpLogCur.write(CRT.join(inputEgIdNotInMgiList))
        fpLogCur.write(CRT + 'Total: %s' % len(inputEgIdNotInMgiList))
    if len(egIdMultiGenesList):
        fpLogCur.write(CRT + CRT + str.center(
            'EG IDs associated with > 1 marker in the Database',60) + CRT)
        fpLogCur.write('%-12s  %-20s  %-20s%s' %
             ('EG ID','PM ID', 'Markers', CRT))
        fpLogCur.write(CRT.join(egIdMultiGenesList))
        fpLogCur.write(CRT + 'Total: %s' % len(egIdMultiGenesList))
    if len(dOrRStatusList):
        fpLogCur.write(CRT + CRT + str.center(
            'References with isDiscard=true or Status=Rejected',60) + CRT)
        fpLogCur.write('%-12s  %-20s  %-20s%s' %
             ('Reference ID','isDiscard', 'Status', CRT))
        fpLogCur.write(CRT.join(dOrRStatusList))
        fpLogCur.write(CRT + 'Total: %s' % len(dOrRStatusList))
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

#
# Main
#

print('%s' % mgi_utils.date())
print('running init')
if init() != 0:
    print('Initialization failed')
    closeFiles()
    sys.exit(1)

print('%s' % mgi_utils.date())
print('running createBCP')
if createBCP() != 0:
    print('Creating BCP Files failed')
    closeFiles()
    sys.exit(1)

print('%s' % mgi_utils.date())
print('running bcpFiles')
if bcpFiles() != 0:
    print('BCP failed')
    closeFiles()
    sys.exit(1)

print('%s' % mgi_utils.date())
print('running updateGoStatus')
if updateGoStatus()  != 0:
    print('Status updates failed')
    closeFiles()
    sys.exit(1)

print('%s' % mgi_utils.date())
print('running writeCuratorLog')
if writeCuratorLog() != 0:
    print('Writing to curator log failed')
    closeFiles()
    sys.exit(1)

closeFiles()

print('%s' % mgi_utils.date())
sys.exit(0)
