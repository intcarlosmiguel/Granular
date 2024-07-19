seed=387
for ((i = 1; i <= 10; i++)); do
    ./main $i 6 60 $seed 0 &
    ((seed++))
    if (($i % 6 == 0)); then
        wait
    fi
done