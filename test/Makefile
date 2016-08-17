CC = h5cc
CFLAGS = -I.. $(DBG)

DBG = -g -DDEBUG

TRASH = *.o *~

all : test11 test12 test13

test11: test11.o ../poh5.o
test12: test12.o ../poh5.o
test13: test13.o ../poh5.o

clean:
	$(RM) $(TRASH)
