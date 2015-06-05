#!/bin/bash

# DATE=$(date +"%Y-%m-%d")
# TIME=$(date +"%T")
# OUTDIR=out_${DATE}_${TIME}
# mkdir $OUTDIR

sed "s|Choice= 1|Choice= 1|" ./dummy.conf > ./config.conf
./run_single
# mv ./time.txt $OUTDIR/time1.txt
# mv ./wealth.txt $OUTDIR/wealth1.txt
# mv ./cooperation.txt $OUTDIR/cooperation1.txt
mv ./time.txt time1.txt
mv ./wealth.txt wealth1.txt
mv ./cooperation.txt cooperation1.txt



sed "s|Choice= 1|Choice= 2|" ./dummy.conf > ./config.conf
./run_single
mv ./time.txt time2.txt
mv ./wealth.txt wealth2.txt
mv ./cooperation.txt cooperation2.txt

sed "s|Choice= 1|Choice= 3|" ./dummy.conf > ./config.conf
./run_single
mv ./time.txt time3.txt
mv ./wealth.txt wealth3.txt
mv ./cooperation.txt cooperation3.txt


sed "s|Choice= 1|Choice= 4|" ./dummy.conf > ./config.conf
./run_single
mv ./time.txt time4.txt
mv ./wealth.txt wealth4.txt
mv ./cooperation.txt cooperation4.txt

# rm config.conf
# mv parameters.txt $OUTDIR/
