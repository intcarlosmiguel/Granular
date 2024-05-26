gcc main.c bib/particle.c bib/vector.c bib/simulate.c bib/mtwister.c -o main -lm -O3
seed=3
for ((i = 1; i <= 500; i++)); do
    # Exemplo de ação: imprimir os valores das variáveis
    ./main 3 3 0 $seed &
    ((seed++))
    if (($i % 3 == 0)); then
        wait
    fi
done

for ((i = 1; i <= 500; i++)); do
    # Exemplo de ação: imprimir os valores das variáveis
    ./main 8 8 0 $seed &
    ((seed++))
    if (($i % 3 == 0)); then
        wait
    fi
done
#python3 plot.py