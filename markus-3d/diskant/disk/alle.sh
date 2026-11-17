set -e
STL=1 python3 diskant.py BESASDES Diskantteil_Markus_v2.step 2>&1 | grep -vE "Warn|warn\("
xvfb-run -a python3 render_d.py d2
for v in ADGC "ADGC - Risch" GCFB GISCISFISH; do n=$(echo "$v" | tr -d ' '); python3 diskant.py "$v" "Diskantteil_Markus_v2_$n.step" 2>&1 | grep -E "reihe|Teile|Error|Trace"; done
