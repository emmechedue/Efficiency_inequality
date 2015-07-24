#!/bin/bash

# DATE=$(date +"%Y-%m-%d")
# TIME=$(date +"%T")
# OUTDIR=out_${DATE}_${TIME}
# mkdir $OUTDIR

sed "s|Choice= 1|Choice= 1|" ./dummy.conf > ./config.conf
./run_ensemble
# mv ./time.txt $OUTDIR/time1.txt
# mv ./wealth.txt $OUTDIR/wealth1.txt
# mv ./cooperation.txt $OUTDIR/cooperation1.txt
mv ./talent.txt talent1.txt
mv ./wealth.txt wealth1.txt
mv ./cooperation.txt cooperation1.txt
mv ./efficiency.txt efficiency1.txt
mv ./giniCoef.txt giniCoef1.txt



sed "s|Choice= 1|Choice= 2|" ./dummy.conf > ./config.conf
./run_ensemble
mv ./talent.txt talent2.txt
mv ./wealth.txt wealth2.txt
mv ./cooperation.txt cooperation2.txt
mv ./efficiency.txt efficiency2.txt
mv ./giniCoef.txt giniCoef2.txt

sed "s|Choice= 1|Choice= 3|" ./dummy.conf > ./config.conf
./run_ensemble
mv ./talent.txt talent3.txt
mv ./wealth.txt wealth3.txt
mv ./cooperation.txt cooperation3.txt
mv ./efficiency.txt efficiency3.txt
mv ./giniCoef.txt giniCoef3.txt


sed "s|Choice= 1|Choice= 4|" ./dummy.conf > ./config.conf
./run_ensemble
mv ./talent.txt talent4.txt
mv ./wealth.txt wealth4.txt
mv ./cooperation.txt cooperation4.txt
mv ./efficiency.txt efficiency4.txt
mv ./giniCoef.txt giniCoef4.txt

# rm config.conf
# mv parameters.txt $OUTDIR/
