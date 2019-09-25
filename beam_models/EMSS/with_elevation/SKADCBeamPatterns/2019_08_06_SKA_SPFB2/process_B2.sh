#!/bin/bash

rm *.mat *.png *.fits interpolated/*.fits
rm -r 2019_08_06_SKA_15_B2
rm -r 2019_08_06_SKA_45_B2
rm -r 2019_08_06_SKA_90_B2

ls -lt

unzip 2019_08_06_SKA_15_B2.zip
unzip 2019_08_06_SKA_45_B2.zip
unzip 2019_08_06_SKA_90_B2.zip

python import_beams_B2.py
mkdir -p interpolated
cd interpolated || exit
python interpolate_beam_B2.py
cd ..
