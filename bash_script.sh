#!/bin/bash
conda info -e

python /run_xgam/Xgam/bin/mkdataselection.py "--config=/run_xgam/Xgam/config/config_dataselection.py"
