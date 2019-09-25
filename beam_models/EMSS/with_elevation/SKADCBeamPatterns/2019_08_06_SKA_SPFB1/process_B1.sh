#!/bin/bash

rm *.mat *.png *.fits interpolated/*.fits
rm -r 2019_08_06_SKA_15_B1
rm -r 2019_08_06_SKA_45_B1
rm -r 2019_08_06_SKA_90_B1

ls -lt

unzip 2019_08_06_SKA_15_B1.zip
unzip 2019_08_06_SKA_45_B1.zip
unzip 2019_08_06_SKA_90_B1.zip

python import_beams_B1.py
mkdir -p interpolated
cd interpolated || exit
python interpolate_beam_B1.py
cd ..
