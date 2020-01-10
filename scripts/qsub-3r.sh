#!/bin/sh

###---Function---sumbit 3r.sh---###

# add all sample names: one conditions for example.
SNAME=(Endosperm_20DAP)

# three conditions will be :
# SNAME=(Endosperm_16DAP Endosperm_20DAP Endosperm_24DAP)

## cycle all of conditions:  i corresponds to numbers of condition. For example "i<10" for 10 conditions
for ((i=0;i<1;i++))
do
   qsub -v SN=${SNAME[$i]} /ingens/home/mintu/scripts/3r.sh
done



