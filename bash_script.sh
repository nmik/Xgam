#!/bin/bash
conda info -e

#python /run_xgam/Xgam/bin/mkdataselection.py "--config=/archive/home/sammazza/fermi_data/CONFIG/config_dataselection_local.py"
python /run_xgam/Xgam/bin/mkmask.py --config=/archive/home/sammazza/fermi_data/CONFIG/config_mask_local.py --srcmask=True --gpmask=flat
