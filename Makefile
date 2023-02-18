a.out: main.c
	gcc main.c

run: a.out
	./a.out

.PHONY: run