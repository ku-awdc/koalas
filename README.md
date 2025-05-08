# koalas
a disease model for Chlamydia in koalas

## TODO

- [X] Replicate results from original paper (deterministic)
- [X] Implement stochastic version
- [] Verify core assumptions for new purpose
- [X] Add compartments Pt, Pf and Pc for test/removed true positive, false positives, and clinically diseased (taking into account se/sp)
- [X] Add active capture for testing/removal strategy (‘Actively’ capturing individuals (at various rate of capture, e.g., 5-6 koalas/day/team); screening by nucleic acid detection (unlikely perfect diagnostic sensitivity and specificity, turn-around time 48-72h); and removing (relocating) positive individuals)
- [X] Add passive capture for testing/remval of clinically diseased (‘Passively’ accessing reported diseased individuals (at various rate of reporting); screening by nucleic acid detection (unlikely perfect diagnostic sensitivity and specificity, turn-around time 48-72h); and removing (relocating) positive individuals)
- [X] Allow natural immunity and vaccination to have different rates - separate V and R
- [] Migrate to C++ for running (not necessary but would be nice)

## Questions / assumptions to assess

- [] Are exponential waiting times reasonable for I and D?
- [] Do we not need an E state?
- [] Vaccination from I/D reverts to susceptible - should this be reversion to I/D?
- [] Which stable baseline parameter set(s) to use as a reference point?
- [] For testing, do we assume that the koala remains in the population (and therefore infectious) while waiting for test results, and then removed after 48-72 hours?
- [] Is it sufficient for vaccine efficacy to be combined with the vaccination rate (i.e. non-successful vaccines simply have no effect?) or do we need to separate efficacy against being infected and efficacy against shedding?
- [] Is it sufficient to assume a single well-mixed group?
- [] Future work:  ABM
