gcc main.c bib/particle.c bib/vector.c bib/simulate.c -o main -lm -O3
time ./main
python3 plot.py