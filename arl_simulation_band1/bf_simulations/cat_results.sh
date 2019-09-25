#!/bin/bash
#!

cat */*_0.csv | head -1 > pointing_simulation_results.csv
cat */*_0.csv */*/*_0.csv | grep -v context >> pointing_simulation_results.csv
