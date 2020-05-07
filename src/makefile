include $(HOME)/make.inc


all: main main.c

main: main.c
	$(CC) $(CFLAGS) $(CVODEINC) $(CAIROINC) -c main.c
	$(CC) -g $(CFLAGS) -o main main.o $(CAIRO) $(CVODE)
	./main

clean:
	rm -f *.o main
