#!/bin/bash

NOWR=$(date --date="yesterday " +"%y%m%d")
cat /u/xavier/idl/xidl/IMG/automated/Reduction/cron/automated_ssh.sh | ssh zdhughes@vhe6.ucsc.edu /bin/bash

mkdir /c/Blazars/data/testbed/nights/$NOWR
cp /usr/local/ftp/incoming/zdhughes/'FLWO_'$NOWR'.tar.gz' /c/Blazars/data/testbed/nights/$NOWR/


