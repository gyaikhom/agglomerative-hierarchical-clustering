agglomerate: agglomerate.c
	gcc -std=c99 -o agglomerate agglomerate.c -lm

clean:
	rm -f agglomerate
