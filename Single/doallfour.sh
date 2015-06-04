#!/bin/bash

DATE=$(date +"%Y-%m-%d")
TIME=$(date +"%T")
OUTDIR=out_${DATE}_${TIME}
mkdir $OUTDIR
PROG=single_run

sed "s|Choice= 1|Choice= 1|" ./dummy.conf > ./config.conf
./$(PROG)
mv ./time.txt $OUTDIR/time1.txt
mv ./wealth.txt $OUTDIR/wealth1.txt
mv ./cooperation.txt $OUTDIR/cooperation1.txt



sed "s|Choice= 1|Choice= 2|" ./dummy.conf > ./config.conf
./$(PROG)
mv ./time.txt $OUTDIR/time2.txt
mv ./wealth.txt $OUTDIR/wealth2.txt
mv ./cooperation.txt $OUTDIR/cooperation2.txt


sed "s|Choice= 1|Choice= 3|" ./dummy.conf > ./config.conf
./$(PROG)
mv ./time.txt $OUTDIR/time3.txt
mv ./wealth.txt $OUTDIR/wealth3.txt
mv ./cooperation.txt $OUTDIR/cooperation3.txt


sed "s|Choice= 1|Choice= 4|" ./dummy.conf > ./config.conf
./$(PROG)
mv ./time.txt $OUTDIR/time4.txt
mv ./wealth.txt $OUTDIR/wealth4.txt
mv ./cooperation.txt $OUTDIR/cooperation4.txt

rm config.conf
mv parameters.txt $OUTDIR/
