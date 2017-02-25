#!/bin/sh
#########################################################################
# sdss_runparallelsrch.sh
# Author: Kathy COoksey            Date: 8 Jun 2011
# Project: SDSS DR7 CIV, SiIV, & CaII search 
# Description: This script parses the large SDSS DR7 QSO structure into
#              multiple parts to divide the labor between processors.
#              The trick is spawn jobs in background and have script
#              just wait for all spawned jobs to finish
#              before continuing.
#
# Call:
#   sdss_runparallelsrch.sh processors list_fil sdsstab
#            sdss-fndlin-options 
#   OR sdss_runparallelsrch.sh processors list_fil sdsstab 
#            cand_fil sdss-civsearch-options] 
# Input:
#   nprocessors -- number of processors; will always run on 1
#   list_fil --  List of QSO spectra to analyze
#   sdsstab -- summary file of QSO properties, in order as list_fil
#   sdss_fndlin-options -- options to be passed to sdss_fndlin
#   cand_fil -- candidate doublets file, if want sdss_civsearch run in
#               parallel.  
#   sdss_civsearch-options -- options to be passed to sdss_civsearch
#   
# Output: 
#   list_fil.#.log -- output from each processor's IDL session from
#                     sdss_fndlin
#   sdss_fndlin's output -- depends on input options
#   cand_fil -- (see above)
#   cand_fil.#.log -- output from each processor's IDL session from
#                     sdss_civsearch
#   cand_fil.cat.log -- output from sdss_catparalleljob
#   These will be stored in whereever the list_fil or cand_fil are
#   stored.
#
# Example: 
# To run sdss_fndlin.pro:
# xterm>sdss_runparallelsrch.sh 4 $XIDL_DIR/SDSS/CIV/dr7qso_srt.list 
#         $XIDL_DIR/SDSS/CIV/dr7qso_srt.fit "/eig,/clobber,pixlim=7,/silent"
# To run sdss_civsearch.pro:
# xterm>sdss_runparallelsrch.sh 4 $XIDL_DIR/SDSS/CIV/dr7qso_srt.list 
#         $XIDL_DIR/SDSS/CIV/dr7qso_srt.fit 
#         $SDSSPATH/$SDSSDR/sdss_siivcand_full.fits "/tcpu,dblt_name='SiIV'"
#
# History:
#   08 Jun 2011  created by KLC
#   09 Jun 2011  Made scalable to any number of processors and have
#                either/or calls
#########################################################################

## Check syntax
if [ $# -lt 4 ] 
then 
    echo "Syntax - sdss_runparallelsrch.sh nprocessors list_file sdsstab "
    echo "               sdss_fndlin-options-required (4 params, for sdss_fndlin)"
    echo "    or - cand_fil sdss_civsearch-options-required (5 params, for sdss_civsearch)"
    exit 1 # fail
fi


if [ $# -eq 4 ]
then
    ######################################################
    ## Run sdss_fndlin.pro: Monitor number of IDL sessions running
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
	echo "Processor $count: sdss_fndlin info piped to $ofil"
	if [ -a $ofil ]
	then rm $ofil
	fi
        ## In order to run an EOF block in the background, must put it in
        ## parentheses with & at the end. Final EOF must be on a line by
        ## itself with no space. Time command prints out run-time
        ## information. 
	(time idl &> $ofil <<EOF
 .r sdss_functions
 sdss_fndlin,'$2','$3',processor=[$count,$1],$4
EOF
	) &
	count=`expr $count + 1`
    done
    
    
    #########
    ## WAIT
    ## Don't run rest of script until IDL sdss_fndlin sessions finish
    wait

    echo ""
    echo "sdss_runparallelsrch.sh: Finished sdss_fndlin call(s)." 
    echo ""
fi



if [ $# -eq 5 ]
then
    ######################################################
    ## Run sdss_civsearch.pro: Monitor number of IDL sessions running
    ## (occupying processors)
    
    count=1
    while [ $count -le $1 ]
    do
        #########
        ## LOOP
        ## Similar to above.
	ofil="$4.$count.log"
	echo "Processor $count: sdss_civsearch info piped to $ofil"
	if [ -a $ofil ]
	then rm $ofil
	fi
	(time idl &> $ofil <<EOF
 .r sdss_functions
 .r sdss_fndlin
 sdss_civsearch,'$2','$3','$4',processor=[$count,$1],$5
EOF
	) &
	count=`expr $count + 1`
    done
    
    
    #########
    ## WAIT
    ## Don't run rest of script until IDL sdss_civsearch sessions finish
    wait
    
    echo ""
    echo "sdss_runparallelsrch.sh: Finished sdss_civsearch call(s)."
    echo ""

    
    ######################################################
    ## CONCATENATE
    ## Then concatenate all the outputs. Doesn't matter what $count is. 
    ofil="$4.cat.log"
    echo "sdss_catparalleljob info piped to $ofil"
    if [ -a $ofil ]
    then rm $ofil
    fi
    (time idl &> $ofil <<EOF
 .r sdss_functions
 sdss_catparalleljob,'$4', [$count,$1]
EOF
    ) &
    
    #########
    ## WAIT
    ## Don't run rest of script until IDL sdss_civsearch sessions finish
    wait

    echo ""
    echo "sdss_runparallelsrch.sh: Finished sdss_catparalleljob call."
    echo ""

fi
