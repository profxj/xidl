#!/bin/sh
#########################################################################
# hyb_runparallelfit.sh
# Author: Kathy COoksey            Date: 8 Jun 2011
# Project: SDSS DR7 CIV, SiIV, & CaII search 
# Description: This script parses the large SDSS DR7 QSO structure into
#              multiple parts to divide the labor between processors.
#              The trick is spawn jobs in background and have script
#              just wait for all spawned jobs to finish
#              before continuing.
#
# Call:
#   hyb_runparallelfit.sh processors list_fil sdsstab
#            sdss-fndlin_fitallspl-options 
#   OR hyb_runparallelfit.sh processors list_fil sdsstab 
#            cand_fil sdss-civsearch-options] 
# Input:
#   nprocessors -- number of processors; will always run on 1
#   list_fil --  List of QSO spectra to analyze
#   sdsstab -- summary file of QSO properties, in order as list_fil
#   sdss_dblfitconti-options -- options to be passed to sdss_dblfitconti
#   cand_fil -- candidate doublets file, if want sdss_civsearch run in
#               parallel.  
#   sdss_civsearch-options -- options to be passed to sdss_civsearch
#   
# Output: 
#   list_fil.#.log -- output from each processor's IDL session from
#                     sdss_dblfitconti
#   sdss_dblfitconti's output -- depends on input options
#   cand_fil -- (see above)
#   cand_fil.#.log -- output from each processor's IDL session from
#                     sdss_civsearch
#   cand_fil.cat.log -- output from sdss_catparalleljob
#   These will be stored in whereever the list_fil or cand_fil are
#   stored.
#
# Example: 
# To run sdss_dblfitconti.pro:
# xterm>hyb_runparallelfit.sh 4 $XIDL_DIR/SDSS/CIV/dr7qso_srt.list 
#         $XIDL_DIR/SDSS/CIV/dr7qso_srt.fit "/eig,/clobber,pixlim=7,/silent"
# To run sdss_civsearch.pro:
# xterm>hyb_runparallelfit.sh 4 $XIDL_DIR/SDSS/CIV/dr7qso_srt.list 
#         $XIDL_DIR/SDSS/CIV/dr7qso_srt.fit 
#         $SDSSPATH/$SDSSDR/sdss_siivcand_full.fits "/tcpu,dblt_name='SiIV'"
#         kcooksey@hawaii.edu
#
# History:
#   08 Jun 2011  created by KLC
#   09 Jun 2011  Made scalable to any number of processors and have
#                either/or calls
#   29 Jul 2011  Adopted from spl_runparallelfit.sh
#   21 Aug 2015  Enable emailed message when finished, KLC
#########################################################################

## Check syntax
if [ $# -lt 4 ] 
then 
    echo "Syntax - hyb_runparallelfit.sh nprocessors list_file sdsstab "
    echo "               sdss_dblfitconti-options-required [email-address]"
    echo "               (4 params, for sdss_dblfitconti)"
    exit 1 # fail
fi


######################################################
## Run sdss_dblfitconti.pro: Monitor number of IDL sessions running
## (occupying processors)
count=1

while [ $count -le $1 ]
do
    #########
    ## LOOP
    ## Start jobs on separate processors in background (assumes new IDL
    ## sessions are farmed to different processors).  Pipe all
    ## information (stdout & stder) from screen to $ofil with &> but must
    ## delete before hand.
    ofil="$2.$count.log.hyb"
    echo "Processor $count: sdss_dblfitconti info piped to $ofil"
    if [ -a $ofil ]
    then rm $ofil
    fi
    ## In order to run an EOF block in the background, must put it in
    ## parentheses with & at the end. Final EOF must be on a line by
    ## itself with no space. Time command prints out run-time
    ## information. 
    (time idl &> $ofil <<EOF
 .r sdss_functions
 .r sdss_dblfitconti
 sdss_dblfitconti,'$2','$3',processor=[$count,$1],$4
EOF
    ) &
    count=`expr $count + 1`
done


#########
## WAIT
## Don't run rest of script until IDL sdss_dblfitconti sessions finish
wait

msg="hyb_runparallelfit.sh: Finished sdss_dblfitconti call(s)."

if [ $# -eq 5 ]
then
    ## Construct and send an email to desired party
    echo $msg > msg.txt
    echo "sdss_dblfitconti,'$2','$3',processor=[*,$1],$4" >> msg.txt
    mail -s 'hyb_runparallelfit.sh' $5 < msg.txt
    \rm msg.txt
fi    

echo ""
echo $msg
echo ""




