#!/bin/bash

set -x
set -e

nx=2000
ny=2000
dx=1
dy=1

mkdir -p bins
rm -f bins/*

cat >config.hrdata <<EOF
{
  nx: $nx
  ny: $ny
  time_steps: 100
  dx: $dx
  dy: $dy
  time_step: 0.3
  cp: cp.bin
  cs: cs.bin
  rho: rho.bin
	save_rate: 100
	save_path: bins
	initial_pressure: initial_pressure.bin
}
EOF

cat >gen_material.py <<EOF
import numpy as np
np.full(($nx, $ny), 2).astype('<f8').tofile("cp.bin")
np.full(($nx, $ny), 1).astype('<f8').tofile("cs.bin")
np.full(($nx, $ny), 1.5).astype('<f8').tofile("rho.bin")
EOF
python gen_material.py

cat >gen_initial.py <<EOF
import numpy as np
data = np.zeros(($nx, $ny), dtype='<f8')
for i in range($nx):
	for j in range($ny):
		if (i - $nx * 0.5) ** 2 + (j - $ny * 0.5) ** 2 < 10 ** 2:
			data[i,j] = 1
data.astype('<f8').tofile("initial_pressure.bin")
EOF
python gen_initial.py

#/usr/bin/time -f "time: %e" ./gc_elastic_2d config.hrdata
/usr/bin/time -f "time: %e" ./gc_elastic_2d_tiled config.hrdata
