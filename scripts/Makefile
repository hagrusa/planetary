CC = gcc
LDFLAGS = -lm


both: CalcP.o CalcBindingE.o
	$(CC) -o CalcP CalcP.o $(LDFLAGS)
	$(CC) -o CalcBindingE CalcBindingE.o $(LDFLAGS)

clean:
	rm -f CalcP.o
	rm -f CalcBindingE.o
#	rm -r CalcP.dSYM
#	rm -r CalcBindingE.dSYM
superclean:
	rm -f CalcP
	rm -f CalcBindingE