#!/bin/sh
#
# pubmed2geneload.sh
###########################################################################
#
#  Purpose:
# 	This script creates pubmed2gene associations in the db
#
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

CONFIG_LOAD=../pubmed2geneload.config

USAGE='Usage: pubmed2geneload.sh"

#
# Make sure the common configuration file exists and source it. 
#
if [ -f ${COMMON_CONFIG} ]
then
    . ${COMMON_CONFIG}
else
    echo "Missing configuration file: ${COMMON_CONFIG}"
    exit 1
fi

#
# Initialize the log file.
# open LOG in append mode and redirect stdout
#
LOG=${LOG_FILE}
rm -rf ${LOG}
>>${LOG}

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
echo "Creating PubMed to Gene Associations" >> ${LOG_DIAG} 2>&1
${PYTHON} ${PUBMED2GENELOAD}/bin/pubmed2geneload.py >> ${LOG_DIAG} 2>&1
STAT=$?
checkStatus ${STAT} "pubmed2geneload.py"

#
# run postload cleanup and email logs
#
shutDown
