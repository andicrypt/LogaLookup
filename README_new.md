What to bench:
Subroutine (including plookup and logalookup) in different threads


Depend on different computer architecture
Subroutine integration with HyperPlonk (including plookup and logalookup) in different threads

Benchmark subroutines
RAYON_NUM_THREADS=8 cargo bench --package subroutines -- lk plk
.
.
RAYON_NUM_THREADS=16 cargo bench --package subroutines -- lk plk
.
.
RAYON_NUM_THREADS=24 cargo bench --package subroutines -- lk plk
.
.
RAYON_NUM_THREADS=32 cargo bench --package subroutines -- lk plk

Benchmark integration with vanilla
RAYON_NUM_THREADS=8 cargo bench --package hyperplonk -- vanilla
.
.
RAYON_NUM_THREADS=16 cargo bench --package hyperplonk -- vanilla
.
.
RAYON_NUM_THREADS=24 cargo bench --package hyperplonk -- vanilla
.
.
RAYON_NUM_THREADS=32 cargo bench --package hyperplonk -- vanilla

Benchmark integration with jelly
RAYON_NUM_THREADS=8 cargo bench --package hyperplonk -- jelly
.
.
RAYON_NUM_THREADS=16 cargo bench --package hyperplonk -- jelly
.
.
RAYON_NUM_THREADS=20 cargo bench --package hyperplonk -- jelly
.
.
RAYON_NUM_THREADS=32 cargo bench --package hyperplonk -- jelly

