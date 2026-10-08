#!/bin/bash
# HB-1a: IHIS time series (every step from postProcessing/hbState) -> one CSV,
# then the comparison figure and metrics against nexaFoam, hy2Foam and dsmcFoam.
# Run from the IHIS case directory; reference CSVs are expected in ./refs/.
set -e

/bin/python3 plot_ttr_vs_time.py \
    --case . \
    --fields TTR TVib_N2 TVib_O2 \
    --output ihis_case3_noVV.csv

/bin/python3 plot_heatbath.py \
    --title "N\$_2\$ V-T relaxation (vib. heating)" \
    --table IHIS     hB_case1.csv \
    --table nexaFoam refs/nexaFoam_case1_heat.csv \
    --ref   hy2Foam  Ttr  refs/hy2Foam_case1_vibHeating_Ttr.csv \
    --ref   hy2Foam  Tvib refs/hy2Foam_case1_vibHeating_Tvib.csv \
    --ref   dsmcFoam Ttr  refs/dsmcFoam_case1_vibHeating_Ttr.csv \
    --ref   dsmcFoam Tvib refs/dsmcFoam_case1_vibHeating_Tvib.csv \
    --teq 7623.3184714709 \
    --output HB1a_comparison.png \
    --metrics HB1a_metrics.csv
