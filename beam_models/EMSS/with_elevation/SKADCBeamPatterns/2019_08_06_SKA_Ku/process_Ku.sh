#!/bin/bash

rm *.mat *.png *.fits interpolated/*.fits
rm -r 2019_08_06_SKA_15_Ku
rm -r 2019_08_06_SKA_45_Ku
rm -r 2019_08_06_SKA_90_Ku

ls -lt

unzip 2019_08_06_SKA_15_Ku.zip
unzip 2019_08_06_SKA_45_Ku.zip
unzip 2019_08_06_SKA_90_Ku.zip

mv Ku_15_12501_5.mat Ku_15_12501.mat
mv Ku_45_12501_5.mat Ku_45_12501.mat
mv Ku_90_12501_5.mat Ku_90_12501.mat
python import_beams_Ku.py
mkdir -p interpolated
cd interpolated || exit
python interpolate_beam_Ku.py
cd ..
