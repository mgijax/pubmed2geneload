
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
#	 3. BCPREF_FILE -  full path to the bcp file
#	 4. BCPSTATUS_FILE -  full path to the bcp file
#	 5. PG_DBUTILS - for the bcp utility script
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
###########################################################################

import sys
import os
import mgi_utils
import db
import loadlib
       
db.setTrace(True)
DEBUG = 0	 # if 0, not in debug mode

# next available MGI_Reference_Assoc, BIB_Workflow_Status primary key
refAssocKey = None
statusKey = None
mgiTypeKey = 2	            # marker
refAssocTypeKey = 1018      # General
createdByKey = 1571         # pubmed2geneload
loaddate = loadlib.loaddate # today's date

# diagnostic log
fpLogDiag = ''

# bcp stuff
outputDir = os.environ['OUTPUTDIR']
refTable = 'MGI_Reference_Assoc'
statusTable = 'BIB_Workflow_Status'
bcpRefFile =  refTable + '.bcp'
bcpStatusFile =  statusTable + '.bcp'
fpRef = ''
fpStatus = ''

# Total gene/reference associations deleted from db
totalDeleted = 0

# references in db that allow status = Indexed updates
statusToDoList = []
statusDoneList = []

# update existing status iscurrent = 0
updateStatusSQL = ""

#
# Purpose: initialize lookups, next primary key, file descriptors
# Returns: 0
# Assumes: input/output files exist and have been opened
# Effects: writes to the file system
# Throws: Nothing
#
def init():
    global refAssocKey, statusKey
    global fpRef, fpStatus, fpLogDiag
    global statusToDoList

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # diagnostic log
    fpLogDiag = open (os.environ['LOG_DIAG'], 'a')

    # bcp file
    fpRef = open(os.environ['BCPREF_FILE'], 'w')
    fpStatus = open(os.environ['BCPSTATUS_FILE'], 'w')

    results = db.sql(''' select nextval('mgi_reference_assoc_seq') as maxKey ''', 'auto')
    refAssocKey = results[0]['maxKey']

    results = db.sql(''' select nextval('bib_workflow_status_seq') as maxKey ''', 'auto')
    statusKey = results[0]['maxKey']

    #
    # reference/pubmedid/geneid associations that already exist in MGI
    #   
    print('getting reference/pubmedid/geneid associations that already exist in MGI')
    db.sql('''
        select c._refs_key, c.pubmedid, a.accid as geneid
        into temporary table currentAssocs
        from MGI_Reference_Assoc ra, BIB_Citation_Cache c, ACC_Accession a
                where c.pubmedid is not null
                and c.jnumid is not null
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
    # GO status = New, Not Routed, Routed, Chosen
    # Relevance = keep
    # Is Reviewed = 0
    # used to update GO status when refs associated with genes
    #
    print('getting list of refs in database that allow status = Indexed update')
    results = db.sql('''
        select c.mgiid
        from BIB_Citation_Cache c, BIB_Workflow_Status s, BIB_Workflow_Relevance v
        where c.jnumid is not null
        and c.isreviewarticle = 0
        and c._Refs_key = s._Refs_key
        and s._Group_key = 31576666
        and s.isCurrent = 1
        and s._Status_key in (31576669, 31576670, 31576671, 71027551)
        and s._Refs_key = v._Refs_key
        and v._Relevance_key = 70594667
        ''', 'auto')
    for r in results:
        statusToDoList.append(r['mgiid'])

    sys.stdout.flush()
    return 0

#
# Purpose: create reference association bcp file, writes discrepancies to data structures by type
# Returns: 0
# Assumes: input/output files exist and have been opened
# Effects: writes to the file system
# Throws: Nothing
#
def createBCP():
    global refAssocKey, statusKey
    global updateStatusSQL

    results = db.sql('''
        select distinct d.geneid, d.pubmedid, a._Object_key as _marker_key, c._refs_key, c.mgiid, c.jnumid
        from DP_EntrezGene_PubMed d, ACC_Accession a, BIB_Citation_Cache c
        where d.taxid = 10090
        and d.geneid = a.accid
        and a._MGIType_key = 2
        and a._LogicalDB_key = 55
        and a.preferred = 1
        and d.pubmedid = c.pubmedid
        and d.jnum is not null
        and not exists (select 1 from currentAssocs c
                where d.pubmedid = c.pubmedid
                and d.geneid = c.geneid
        )
        and not exists (select 1 from pmMultiples p where d.pubmedid = p.pubmedid)
        and not exists (select 1 from egMultiples p where d.geneid = p.egid)
        ''', 'auto')

    for r in results:

        pmID = r['pubmedid']
        markerKey = r['_marker_key']
        refKey = r['_refs_key']
        refid = r['mgiid']
        jnumid = r['jnumid']

        # if GO/Status rules are met
        #       add to reference/marker assoc
        #       add to GO/Status
        if refid in statusToDoList:
                if refid not in statusDoneList:
                        updateStatusSQL += 'update BIB_Workflow_Status set iscurrent = 0 where _refs_key = %s and _group_key = 31576666 and _status_key in (31576669, 31576670, 31576671, 71027551) and iscurrent = 1;\n' % (refKey)
                        fpStatus.write('%s|%s|31576666|31576673|1|%s|%s|%s|%s\n' % \
                                (statusKey, refKey, createdByKey, createdByKey, loaddate, loaddate))
                        statusKey = statusKey + 1
                        statusDoneList.append(refid)

        fpRef.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' % \
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

    fpRef.close()
    fpStatus.close()
    #print(updateStatusSQL)
    db.commit()

    if DEBUG:
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
        (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), refTable, outputDir, bcpRefFile)
    print(bcpCmd)
    fpLogDiag.write('%s\n' % bcpCmd)
    print('executing bcp')
    os.system(bcpCmd)

    bcpCmd = '%s %s %s %s %s %s "\|" "\\n" mgd' % \
        (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), statusTable, outputDir, bcpStatusFile)
    print(bcpCmd)
    fpLogDiag.write('%s\n' % bcpCmd)
    print('executing bcp')
    os.system(bcpCmd)

    db.sql(updateStatusSQL)
    db.commit()

    # update mgi_reference_assoc auto-sequence
    db.sql(''' select setval('mgi_reference_assoc_seq', (select max(_assoc_key) from MGI_Reference_Assoc)) ''', None)
    db.commit()
    db.sql(''' select setval('bib_workflow_status_seq', (select max(_assoc_key) from BIB_Workflow_Status)) ''', None)
    db.commit()

    return 0

#
# Purpose: Close files.
# Returns: 0
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def closeFiles():
    if fpRef:
        fpRef.close()
    if fpStatus:
        fpStatus.close()
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

closeFiles()
print('%s' % mgi_utils.date())
sys.exit(0)

