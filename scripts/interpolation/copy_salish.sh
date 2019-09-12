#!/bin/bash

# Script to copy a set of files from Susan Allen's SalishSeaCast model

start=2017-04-01
end=2017-05-01

startdate=$(date -I -d $start)
enddate=$(date -I -d $end)

RESULTS=/results/SalishSea/nowcast-green
SAVEDIR=/data/hdd/salishsea/model_orig

d=$startdate
while [ "$d" != "$enddate" ]; do
    echo $d
    datedir_local=$(date -d $d +'%Y%m%d')
    if [[ ! -d $SAVEDIR/$datedir_local ]]; then
       mkdir $SAVEDIR/$datedir_local
    fi
    datedir_salish="$(date -d $d +'%d%b%y')"
    datedir_salish="${datedir_salish,,}"
    echo $datedir_salish
    scp salish:$RESULTS/$datedir_salish/*1h*_grid_T.nc $SAVEDIR/$datedir_local/.
    scp salish:$RESULTS/$datedir_salish/*1h*_grid_W.nc $SAVEDIR/$datedir_local/.
    scp salish:$RESULTS/$datedir_salish/*1h*_grid_U.nc $SAVEDIR/$datedir_local/.
    scp salish:$RESULTS/$datedir_salish/*1h*_grid_V.nc $SAVEDIR/$datedir_local/.
    scp salish:$RESULTS/$datedir_salish/*1h*_ptrc_T.nc $SAVEDIR/$datedir_local/.
    d=$(date -I -d "$d + 1 day")
done
