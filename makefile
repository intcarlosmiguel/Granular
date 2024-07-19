main: main.c
	gcc main.c -o main -lm

clean:
	rm main

run: main
	./main 6 6 60 7152 0