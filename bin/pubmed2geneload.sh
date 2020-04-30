#!/bin/sh
#
# pubmed2geneload.sh
###########################################################################
#
#  Purpose:
# 	This script creates pubmed2gene associations in the db
#
  Usage=pubmed2geneload.sh
#
#  Env Vars:
#
#      See the configuration file
#
#  Inputs:
#
#      - Common configuration file -
#               /usr/local/mgi/live/mgiconfig/master.config.sh
#      - Load configuration file - pubmed2geneload.config
#      - RADAR database
#
#  Outputs:
#
#      - An archive file
#      - Log files defined by the environment variables ${LOG_PROC},
#        ${LOG_DIAG}, ${LOG_CUR} and ${LOG_VAL}
#      - Records written to the database tables
#      - Exceptions written to standard error
#      - Configuration and initialization errors are written to a log file
#        for the shell script
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#      2:  Non-fatal error occurred
#
#  Assumes:  Nothing
#
# History:
#
# sc	09/04/2017 - created
#

cd `dirname $0`
LOG=`pwd`/pubmed2geneload.log
rm -rf ${LOG}

CONFIG_LOAD=../pubmed2geneload.config

#
# verify & source the configuration file
#

if [ ! -r ${CONFIG_LOAD} ]
then
    echo "Cannot read configuration file: ${CONFIG_LOAD}"
    exit 1
fi

. ${CONFIG_LOAD}

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi
echo ${MGD_DBSERVER}
echo ${MGD_DBNAME}

#####################################
#
# Main
#
#####################################
#
# createArchive including OUTPUTDIR, startLog, getConfigEnv
# sets "JOBKEY"

preload ${OUTPUTDIR} 

#
# Create associations in the database
#
echo "Creating PubMed to Gene Associations" >> ${LOG_DIAG}
${PYTHON} ${PUBMED2GENELOAD}/bin/pubmed2geneload.py >> ${LOG_DIAG}
STAT=$?
checkStatus ${STAT} "pubmed2geneload.py"

#
# run postload cleanup and email logs
#
shutDown
