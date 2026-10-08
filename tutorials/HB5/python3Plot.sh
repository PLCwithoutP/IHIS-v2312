#!/bin/bash
/bin/python3 plot_ttr_vs_time.py \
	--case "hB_case5_TTR" .  \
	--field TTR \
	--title "Air Chemical Relaxation" \
	--output air_ttr.png \
	--export-csv

