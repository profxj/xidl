#!/bin/sh
#########################################################################
# eig_runparallelfit.sh
# Author: Kathy COoksey            Date: 8 Jun 2011
# Project: SDSS DR7 CIV, SiIV, & CaII search 
# Description: This script parses the large SDSS DR7 QSO structure into
#              multiple parts to divide the labor between processors.
#              The trick is spawn jobs in background and have script
#              just wait for all spawned jobs to finish
#              before continuing.
#
# Call:
#   eig_runparallelfit.sh processors sdsstab runeigqsoconti-options 
#
# Input:
#   nprocessors -- number of processors; will always run on 1
#   sdsstab -- summary file of QSO properties, in order as list_fil
#   runeigqsoconti-options -- options to be passed to runeigqsoconti
#   
# Output: 
#   list_fil.#.log -- output from each processor's IDL session from
#                     runeigqsoconti
#   runeigqsoconti's output -- depends on input options
#   These will be stored in whereever the list_fil is stored.
#
# Example: 
# To run runeigqsoconti.pro:
# xterm>eig_runparallelfit.sh 4 $XIDL_DIR/SDSS/CIV/dr7qso_srt.fit 
#         "/clobber" kcooksey@hawaii.edu
#
# History:
#   08 Jun 2011  created by KLC
#   09 Jun 2011  Made scalable to any number of processors and have
#                either/or calls
#   27 Jul 2011 Adopted from sdss_runparallelsrch.sh
#   21 Aug 2015  Enable emailed message when finished, KLC
#########################################################################

## Check syntax
if [ $# -lt 3 ] 
then 
    echo "Syntax - eig_runparallelfit.sh nprocessors sdsstab "
    echo "               runeigqsoconti-options-required [email-address]"
    echo "               (3 params, for runeigqsoconti)"
    exit 1 # fail
fi


######################################################
## Run runeigqsoconti.pro: Monitor number of IDL sessions running
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
    ofil="$2.$count.log.eig"
    echo "Processor $count: runeigqsoconti info piped to $ofil"
    if [ -a $ofil ]
    then rm $ofil
    fi
    ## In order to run an EOF block in the background, must put it in
    ## parentheses with & at the end. Final EOF must be on a line by
    ## itself with no space. Time command prints out run-time
    ## information. 
    (time idl &> $ofil <<EOF
 .r sdss_functions
 runeigqsoconti,sdsssum='$2',processor=[$count,$1],$3
EOF
    ) &
    count=`expr $count + 1`
done

    
#########
## WAIT
## Don't run rest of script until IDL runeigqsoconti sessions finish
wait

msg="eig_runparallelfit.sh: Finished runeigqsoconti call(s)."

if [ $# -eq 4 ]
then
    ## Construct and send an email to desired party
    echo $msg > msg.txt
    echo "runeigqsoconti,sdsssum='$2',processor=[*,$1],$3" >> msg.txt
    mail -s 'eig_runparallelfit.sh' $4 < msg.txt
    \rm msg.txt
fi

echo ""
echo $msg
echo ""




