#!/bin/sh
#########################################################################
# sdss_runparallelcmplt.sh
# Author: Kathy COoksey            Date: 8 Jun 2011
# Project: SDSS DR7 CIV, SiIV, & CaII search 
# Description: This script parses the large SDSS DR7 QSO structure into
#              multiple parts to divide the labor between processors.
#              The trick is spawn jobs in background and have script
#              just wait for all spawned jobs to finish
#              before continuing.
#
# Call:
#   sdss_runparallelcmplt.sh processors list_fil sdsstab
#            sdss-completeness-options 
#
# Input:
#   nprocessors -- number of processors; will always run on 1
#   list_fil --  List of QSO spectra to analyze
#   sdsstab -- summary file of QSO properties, in order as list_fil
#   sdss_completeness-options -- options to be passed to sdss_completeness
#   cand_fil -- candidate doublets file, if want sdss_civsearch run in
#               parallel.  
#   sdss_civsearch-options -- options to be passed to sdss_civsearch
#   
# Output: 
#   list_fil.#.log -- output from each processor's IDL session from
#                     sdss_completeness
#   sdss_completeness's output -- depends on input options
#   cand_fil -- (see above)
#   cand_fil.#.log -- output from each processor's IDL session from
#                     sdss_civsearch
#   cand_fil.cat.log -- output from sdss_catparalleljob
#   These will be stored in whereever the list_fil or cand_fil are
#   stored.
#
# Example: 
# To run sdss_completeness.pro:
# xterm>sdss_runparallelcmplt.sh 4 $XIDL_DIR/SDSS/CIV/dr7qso_srt.list 
#         $XIDL_DIR/SDSS/CIV/dr7qso_srt.fit "/eig,/clobber,pixlim=7,/silent"
#         kcooksey@hawaii.edu
#
# History:
#   08 Jun 2011  created by KLC
#   09 Jun 2011  Made scalable to any number of processors and have
#                either/or calls
#   21 Dec 2011  Adopted from sdss_runparallelsrch.sh, KLC
#   21 Aug 2015  Enable emailed message when finished, KLC
#########################################################################

## Check syntax
if [ $# -lt 5 ] 
then 
    echo "Syntax - sdss_runparallelcmplt.sh nprocessors list_file sdsstab nabstot"
    echo "               sdss_completeness-options-required [email-address]"
    echo "               (5 params, for sdss_completeness)"

    exit 1 # fail
fi



######################################################
## Run sdss_completeness.pro: Monitor number of IDL sessions running
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
    ofil="$2.$count.log"
    echo "Processor $count: sdss_completeness info piped to $ofil"
    if [ -a $ofil ]
    then rm $ofil
    fi
    ## In order to run an EOF block in the background, must put it in
    ## parentheses with & at the end. Final EOF must be on a line by
    ## itself with no space. Time command prints out run-time
    ## information. 
    (time idl &> $ofil <<EOF
 .r sdss_functions
 sdss_completeness,'$2','$3','$4',processor=[$count,$1],$5
EOF
    ) &
    count=`expr $count + 1`
done


#########
## WAIT
## Don't run rest of script until IDL sdss_completeness sessions finish
wait

msg="sdss_runparallelcmplt.sh: Finished sdss_completeness call(s)."

if [ $# -eq 6 ] 
then
    ## Construct and send an email to desired party
    echo $msg > msg.txt
    echo "sdss_completeness,'$2','$3','$4',processor=[*,$1],$5" >> msg.txt
    mail -s 'sdss_runparallelcmplt.sh' $6 < msg.txt
    \rm msg.txt
fi

echo ""
echo $msg
echo ""




