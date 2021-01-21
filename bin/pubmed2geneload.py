
##########################################################################
#
# Purpose:
#       Create pubmed to gene associations in the database
#
# Usage: pubmed2geneload.py
#
# Env Vars:
#	 1. OUTPUTDIR - where to write the bcp file
#	 2. LOG_DIAG - diagnostic log 
#	 3. BCP_FILE -  full path to the bcp file
#	 4. PG_DBUTILS - for the bcp utility script
#	 5. JAVA_API_URL - access to the java api for status updates
#	 6. JAVA_API_TOKEN - security token for java api updates
#	 7. UPDATE_BATCH - size of the status update batches
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
#
#       01/20/2021      lec
#       - TR13349/Genome Build 39 project
#
#	2/9/2018	sc
#	- TR12760 - create marker/reference associations and update GO status
#		only if reference has JNumber
###########################################################################

import sys
import os
import mgi_utils
import db
import loadlib
import subprocess
        
db.setTrace(True)

DEBUG = 0	 # if 0, not in debug mode
bcpon = 1        # bcp into the database?  default is yes.

# for creating API URL to update status
PUT = 'curl -X PUT "'
JAVA_API_URL = os.environ['JAVA_API_URL']
JAVA_API_TOKEN = os.environ['JAVA_API_TOKEN']
STATUS_URL = 'littriage/statusUpdate?accid=%s&group=GO&status=Indexed" -H "accept: application/json" -H  "api_access_token: ' + JAVA_API_TOKEN + '" -H  "username: pm2geneload"'

# the full URL - just plug in comma delim (no space) list of reference MGI ID
FULL_API_URL = PUT + JAVA_API_URL + STATUS_URL

# max number of status updates in a batch
UPDATE_BATCH = int(os.environ['UPDATE_BATCH'])

# next available MGI_Reference_Assoc primary key
refAssocKey = None
mgiTypeKey = 2	            # marker
refAssocTypeKey = 1018      # General
createdByKey = 1571         # pubmed2geneload
loaddate = loadlib.loaddate # today's date

# diagnostic log
fpLogDiag = ''

# bcp stuff
outputDir = os.environ['OUTPUTDIR']
refTable = 'MGI_Reference_Assoc'
bcpFile =  refTable + '.bcp'
fpBcp = ''

# Total gene/reference associations deleted from db
totalDeleted = 0

# references in db that allow status = Indexed updates
dbStatusRefList = []

# references from associations that allow status = Indexed updates
assocStatusRefList = []

#
# Purpose: initialize lookups, next primary key, file descriptors
# Returns: 0
# Assumes: input/output files exist and have been opened
# Effects: writes to the file system
# Throws: Nothing
#
def init():
    global refAssocKey, fpBcp, fpLogDiag
    global dbStatusRefList

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # diagnostic log
    fpLogDiag = open (os.environ['LOG_DIAG'], 'a')

    # bcp file
    fpBcp = open(os.environ['BCP_FILE'], 'w')

    results = db.sql(''' select nextval('mgi_reference_assoc_seq') as maxKey ''', 'auto')
    refAssocKey = results[0]['maxKey']

    #
    # reference/pubmedid/geneid associations that already exist in MGI
    #   
    print('getting reference/pubmedid/geneid associations that already exist in MGI')
    db.sql('''
        select c._refs_key, c.pubmedid, a.accid as geneid
        into temporary table currentAssocs
        from MGI_Reference_Assoc ra, BIB_Citation_Cache c, ACC_Accession a
                where c.pubmedid is not null
                and c._refs_key = ra._refs_key
                and ra._MGIType_key = 2
                and ra._refassoctype_key = 1018
                and ra._CreatedBy_key != 1571
                and ra._Object_key = a._Object_key 
                and a._MGIType_key = 2
                and a._LogicalDB_key = 55
                and a.preferred = 1
    ''', None)
    db.sql('create index associdx1 on currentAssocs(_refs_key)', None)
    db.sql('create index associdx2 on currentAssocs(pubmedid)', None)
    db.sql('create index associdx3 on currentAssocs(geneid)', None)

    #
    # pubmed -> eg group > 15
    #
    print('getting pubmed -> eg group > 15')
    db.sql('''
        select pubmedid
        into temporary table pmMultiples
        from DP_EntrezGene_PubMed
        where taxid = 10090
        group by pubmedid having count(*) > 15
        ''', None)
    db.sql('create index pmidx1 on pmMultiples(pubmedid)', None)
        
    #
    # egid with > 1 marker
    #
    print('getting egid with > 1 marker')
    results = db.sql('''
        WITH egid AS (
        select a.accid as egid
        from ACC_Accession a
        where a._MGIType_key = 2
        and a._LogicalDB_key = 55
        and a.preferred = 1
        group by accid having count(*) > 1 
        )
        select egid
        into temporary table egMultiples
        from egid, ACC_Accession a
        where egid = a.accID
        and a._MGIType_key = 2
        and a._LogicalDB_key = 55
        and a.preferred = 1
        ''', 'auto')
    db.sql('create index egidx1 on egMultiples(egid)', None)

    #
    # GO status = Not Routed, Routed, Chosen
    # Relevance = keep
    # Is Reviewed = 0
    # used to update GO status when refs associated with genes
    #
    print('getting list of refs in database that allow status = Indexed update')
    results = db.sql('''
        select c.mgiid
        from BIB_Citation_Cache c, BIB_Workflow_Status s, BIB_Workflow_Relevance v
        where c._Refs_key = s._Refs_key
        and c.isreviewarticle = 0
        and s._Group_key = 31576666
        and s.isCurrent = 1
        and s._Status_key in (71027551, 31576670, 31576671, 71027551)
        and s._Refs_key = v._Refs_key
        and v._Relevance_key = 70594667
        ''', 'auto')
    for r in results:
        dbStatusRefList.append(r['mgiid'])

    sys.stdout.flush()
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
    global refAssocKey
    global assocStatusRefList

    results = db.sql('''
        select distinct d.geneid, d.pubmedid, a._Object_key as _marker_key, c._refs_key, c.mgiid, c.jnumid
        from DP_EntrezGene_PubMed d, ACC_Accession a, BIB_Citation_Cache c
        where d.taxid = 10090
        and d.geneid = a.accid
        and a._MGIType_key = 2
        and a._LogicalDB_key = 55
        and a.preferred = 1
        and d.pubmedid = c.pubmedid
        and not exists (select 1 from currentAssocs c
                where d.pubmedid = c.pubmedid
                and d.geneid = c.geneid
        )
        and not exists (select 1 from pmMultiples p where d.pubmedid = p.pubmedid)
        and not exists (select 1 from egMultiples p where d.geneid = p.egid)
        ''', 'auto')

    for r in results:

        addToBcp = 0

        pmID = r['pubmedid']
        markerKey = r['_marker_key']
        refKey = r['_refs_key']
        refid = r['mgiid']
        jnumid = r['jnumid']

        # if jnum exists, add to reference/marker assoc
        if jnumid != None:
                addToBcp = 1

        # if GO/Status rules are met
        #       add association to reference/marker
        #       add to GO/Status
        if refid in dbStatusRefList:
                addToBcp = 1
                assocStatusRefList.append(refid)

        # if jnum exists, add to reference/marker assoc
        # or
        # if GO/Status rules are met
        if addToBcp == 1:

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
    global totalDeleted

    fpBcp.close()
    db.commit()

    if DEBUG or not bcpon:
        return 0

    # delete from MGI_Reference_Assoc
    print('deleting from MGI_Reference_Assoc')
    db.sql('''
        select _Assoc_key
        into temporary table toDelete
        from MGI_Reference_Assoc
        where _CreatedBy_key = 1571
        ''', None)

    db.sql('''create index idx1 on toDelete(_Assoc_key)''', None)

    results = db.sql('''select count(*) as deleteCt
        from toDelete''', 'auto')

    for r in results:
        totalDeleted = r['deleteCt']

    db.sql('''
        delete from MGI_Reference_Assoc a
        using toDelete d
        where a._Assoc_key = d._Assoc_key
        ''', None)

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
# Assumes:  assocStatusRefList has been loaded
# Effects:  Updates GO status in the database
# Throws:   nothing
#
def updateGoStatus():
    global assocStatusRefList

    # for all reference associations we are creating:
    #       change status = Indexed
    #       java api will generate jnum, if necessary

    # get list of refids that actually need updating, reporting as we go
    print('updateGoStatus size of assocStatusRefList: %s' % len(assocStatusRefList))

    while assocStatusRefList != []:

        batchToRun = ','.join(assocStatusRefList[0:UPDATE_BATCH])
        print(batchToRun)
        checkGoStatus(batchToRun)
        del assocStatusRefList[0:UPDATE_BATCH]

        print('%s' % mgi_utils.date())
        print('running subprocess')
        print(FULL_API_URL % batchToRun)

        cmd = FULL_API_URL % batchToRun
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        stdout = result.stdout
        stderr = result.stderr
        returnCode = result.returncode

        print('after subprocess stdout: %s stderr: %s returnCode: %s' % (stdout, stderr, returnCode))
        checkGoStatus(batchToRun)

        if returnCode != 0:
            return 1

    return 0

def checkGoStatus(batchToRun):

    testbatchToRun = batchToRun.replace(",MGI:", "',MGI:")
    testbatchToRun = testbatchToRun.replace("MGI:", "'MGI:")
    testbatchToRun += "'"

    results = db.sql('''
        select c._refs_key, c.mgiID, c.jnumid, c.pubmedid, t.term as relvance, t2.term as status,
        u.login as relevanceuser, u2.login as littriageuser, c.short_citation
        from bib_citation_cache c, bib_refs r, bib_workflow_relevance v, voc_term t, mgi_user u,
                bib_workflow_status s, voc_term t2, voc_term t3, mgi_user u2
        where r._refs_key = c._refs_key
        and r._refs_key = v._refs_key
        and v.isCurrent = 1
        and v._relevance_key = 70594667
        and v._relevance_key = t._term_key
        and v._createdby_key = u._user_key
        and r._refs_key = s._refs_key
        and s.isCurrent = 1
        and s._Group_key = 31576666
        and s._status_key = t2._term_key
        and s._group_key = t3._term_key
        and s._createdby_key = u2._user_key
        and c.mgiid in (%s)
        order by c.short_citation, t3.term
        ''' % (testbatchToRun), 'auto')
    for r in results:
        print(r)

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
    if fpLogDiag:
        fpLogDiag.close()
    return 0

#
# Main
#

print('%s' % mgi_utils.date())
print('running init')
sys.stdout.flush()
if init() != 0:
    print('Initialization failed')
    closeFiles()
    sys.exit(1)

print('%s' % mgi_utils.date())
print('running createBCP')
sys.stdout.flush()
if createBCP() != 0:
    print('Creating BCP Files failed')
    closeFiles()
    sys.exit(1)

print('%s' % mgi_utils.date())
print('running bcpFiles')
sys.stdout.flush()
if bcpFiles() != 0:
    print('BCP failed')
    closeFiles()
    sys.exit(1)

print('%s' % mgi_utils.date())
print('running updateGoStatus')
sys.stdout.flush()
if updateGoStatus()  != 0:
    print('Status updates failed')
    closeFiles()
    sys.exit(1)

closeFiles()

print('%s' % mgi_utils.date())
sys.exit(0)
