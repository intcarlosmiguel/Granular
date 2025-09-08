gcc main.c -o main -lm -O3 -fopenmp
seed=38
for abertura in $(seq 45 5 60); do
    ./main 1 1 60 $seed 0 50 $abertura
    seed=$((seed + 500))
done