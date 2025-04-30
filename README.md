# koalas
a disease model for Chlamydia in koalas

## TODO

[] Replicate results from original paper (deterministic)
[] Implement stochastic version
[] Verify core assumptions for new purpose
[] Add active capture for testing/removal strategy (‘Actively’ capturing individuals (at various rate of capture, e.g., 5-6 koalas/day/team); screening by nucleic acid detection (unlikely perfect diagnostic sensitivity and specificity, turn-around time 48-72h); and removing (relocating) positive individuals)
[] Add passive capture for testing/remval of clinically diseased (‘Passively’ accessing reported diseased individuals (at various rate of reporting); screening by nucleic acid detection (unlikely perfect diagnostic sensitivity and specificity, turn-around time 48-72h); and removing (relocating) positive individuals)
[] Allow natural immunity and vaccination to have different rates
[] Migrate to C++ for running (not necessary but would be nice)

## Questions / assumptions to assess

[] Is it reasonable to vaccinate I/D animals?
[] Treatment from I/D reverts to susceptible, not R?
[] Which stable baseline parameter set(s) to use as a reference point?
[] Is it sufficient for vaccine efficacy to be combined with the vaccination rate (i.e. non-successful vaccines simply have no effect?) or do we need to separate efficacy against being infected and efficacy against shedding?
[] Is it sufficient to assume a single well-mixed group?
