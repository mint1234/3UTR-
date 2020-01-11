#!/bin/sh

###---Function---sumbit 3r.sh---###

# add all sample names: two conditions for example.
SNAME=(Endosperm_16DAP Endosperm_20DAP)

## cycle all of conditions:  i corresponds to numbers of condition. 
for ((i=0;i<2;i++))
do
   qsub -v SN=${SNAME[$i]} /ingens/home/mintu/scripts/3r.sh
done



