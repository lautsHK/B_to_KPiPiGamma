#!/bin/bash
charge="1" # Only 0 or 1
sqrt_s_min="0.778" # At least larger than sum of the masses (pi + pi + K), i.e. 0.778 GeV 
sqrt_s_max="1.6" # At least larger than sum of the masses (pi + pi + K), i.e. 0.778 GeV, and please larger than sqrt_s_min
mode="1" # 0 for constant sqrt_s, 1 for sqrt_s_min < sqrt_s < sqrt_s_max.
#ss1270="1.36603"
ss1270="1.36603" # signal strength of K1270
#ss1400="-0.366025"
ss1400="-0.366025"
ss1410="0."
ss1430="0."
ss1680="0."
lambda="1.0" # lambda should be between -1.0 and 1.0 (1.0 for fully right hand polarised, and -1.0 for left hand.)
n_events="500"
jobname="FAST"

echo "Running job $jobname with parameters:"
echo "  charge: $charge"
echo "  sqrt_s_min: $sqrt_s_min"
echo "  sqrt_s_max: $sqrt_s_max"
echo "  mode: $mode"
echo "  ss1270: $ss1270"
echo "  ss1400: $ss1400"
echo "  ss1410: $ss1410"
echo "  ss1430: $ss1430"
echo "  ss1680: $ss1680"
echo "  lambda: $lambda"
echo "  n_events: $n_events"

cd ./../build
./run_B_to_KPiPiGamma_FAST "$charge" "$sqrt_s_min" "$sqrt_s_max" "$mode" "$ss1270" "$ss1400" "$ss1410" "$ss1430" "$ss1680" "$lambda" "$n_events" "$jobname"

