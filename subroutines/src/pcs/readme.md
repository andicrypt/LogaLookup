KZG based multilinear polynomial commitment
-----

### Compiling features:
- `parallel`: use multi-threading when possible.
- `print-trace`: print out user friendly information about the running time for each micro component.
- `extensive_sanity_checks`: runs additional sanity checks that is not essential and will slow down the scheme.


### Running Benchmark for LogaLookup and Plookup
```
    RAYON_NUM_THREADS=<number_of_threads> cargo bench --package subroutines -- lk plk
```

### Running Benchmark for integration of LogaLookup or Plookup into HyperPlonk
* For Jelly gate:
```
    RAYON_NUM_THREADS=<number_of_threads> cargo bench --package hyperplonk -- jelly
```
* For Vanilla gate:
```
    RAYON_NUM_THREADS=<number_of_threads> cargo bench --package hyperplonk -- vanilla
```

### Drawing comparison chart for Benchmarking
```
python3 plotter.py
```
