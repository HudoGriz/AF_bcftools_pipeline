#!/bin/bash

# Replace all "${log_file}" with "${base_name}.log" in the script sections
sed -i 's/"\${log_file}"/"\${base_name}.log"/g' main.nf

# Replace val(log_file) with path(log) in inputs/outputs
sed -i 's/val(log_file)/path(log)/g' main.nf
sed -i 's/, val(log_file)/, path(log)/g' main.nf

echo "Fixed all log_file references"
