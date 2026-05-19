# CACTUS TREFMAX Greedy Sweep

This directory is the final omnetpp_17, TRHD=250, rate8 sweep for the
REFab-slack greedy CACTUS mode.

Command:

```bash
JOBS=8 ./run_trefmax_sweep.sh
```

Key files:

- `trefmax_sweep_omnetpp_trhd250.csv`: compact summary for direct/random and
  `tref_interval = 1, 2, 4, 8`.
- `mop4_cactus_fig15_*_trefmax*_trhd250_omnetpp_17.mtf.gz.out`: raw simulator
  outputs used by the CSV.

Current takeaway:

- `random` with interval 1/2/4 reaches avg IPC 0.753 and eliminates
  ABO-triggered RFMab for this trace.
- `direct` needs interval 1/2 to eliminate almost all alerts/RFMab; interval
  4/8 lets RFMab return and IPC drops.
- The run script removes NUL bytes from temporary stdout before publishing each
  `.out`, because parallel progress printing can otherwise leave terminal NUL
  artifacts in the text output.
