#!/bin/bash

python /run_xgam/Xgam/bin/mkdataselection.py --config=/run_xgam/home/CONFIG/config_dataselection.py

python /run_xgam/Xgam/bin/mkforeground.py --infile=/run_xgam/Xgam/fits/gll_iem_v07.fits --nsideout=1024

python /run_xgam/Xgam/bin/mksmartmask.py --config=/run_xgam/home/CONFIG/config_datafluxmaps.py --srccat=/run_xgam/Xgam/fits/gll_psc_v20.fit --srcextcat=/run_xgam/Xgam/fits/gll_psc_v20.fit --irfs=P8R3_SOURCEVETO_V2 --evtype=3 --ltfile=/run_xgam/home/output/output_gtltcube/tag_outofltcube.fits --outflabel=tag

python /run_xgam/Xgam/bin/mkdatafluxmaps.py --config=/run_xgam/home/CONFIG/config_datafluxmaps.py
