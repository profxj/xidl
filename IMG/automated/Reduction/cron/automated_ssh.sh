#!/bin/bash

NOW=$(date --date="yesterday " +"20%y.%m%d")
NOWR=$(date --date="yesterday " +"%y%m%d")
mkdir /home/vhep/zdhughes/tmp/$NOWR
cp /data_disks/vhe6/48inch_data/2012/$NOW/* /home/vhep/zdhughes/tmp/$NOWR
cd /home/vhep/zdhughes/tmp/$NOWR
tar -zcvf 'FLWO_'$NOWR'.tar.gz' *fits*
#rm *fits*
rm *ccd*


HOST=ftp.ucolick.org  #This is the FTP servers host or IP address.
USER=anonymous             #This is the FTP user that has access to the server.
PASS=zdhughes@ucolick.org          #This is the password for the FTP user.

# Call 1. Uses the ftp command with the -inv switches.  -i turns off interactive prompting. -n Restrains FTP from attempting the auto-login feature. -v enables verbose and progress. 

ftp -inv $HOST << EOF

# Call 2. Here the login credentials are supplied by calling the variables.

user $USER $PASS

# Call 3. Here you will change to the directory where you want to put or get
cd /pub/incoming/zdhughes/

# Call4.  Here you will tell FTP to put or get the file.
put *tar*

# or
#get test.txt
bye


EOF  
