#!/bin/bash

##
## @file run_TMPL.sh
##
## Batch script template for mogon, process by your script
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2014-12-02

# Setup of BSUB
#BSUB -J MailJobIsDone4_DYTURBO
#BSUB -L /bin/bash
#BSUB -B
#BSUM -u cuth@uni-mainz.de
#BSUB -o /etapfs03/atlashpc/cuth/mail.log
#BSUB -e /etapfs03/atlashpc/cuth/mail.err
#BSUB -q atlasshort
#BSUB -W 3:00
#BSUB -app Reserve500M
#BSUB -n 1
#BSUB -R 'rusage[atlasio=0]'
#BSUB -w 'ended(dyturbo_*)' 

# ===================================
# just send email that everything is finished
text="Hey listna,\n all jobs \`dyturbo_*\` have been done at `date`.\n\nYours,\nMogon"
echo -e $text | mail -s "Mogon job report" cuth@uni-mainz.de
exit 0


