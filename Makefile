CC=gcc

#flags de compilacao
CFLAGS= -g -Wall -lm -I/usr/local/include -L/usr/local/lib -lmatheval
#arquivos-objeto
OBJ = newtonSNL.o

#regras de compilacao
newtonSNL = newtonSNL.o utils.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) 

#regra default
all: newtonSNL

#regras de ligacao
newtonSNL: $(newtonSNL)

	$(CC) -o $@ $^ $(CFLAGS) 

.PHONY:	clean

#remove arquivos temporarios
clean: 

	rm *.o 
