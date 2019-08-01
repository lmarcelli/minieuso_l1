#################################################################
####   Makefile to compile Level1 Mini-EUSO Software   ##########
#################################################################
euso_l1: euso_l1.o
	gcc -o euso_l1 euso_l1.o `root-config --libs --cflags` -I `root-config --incdir`  -lstdc++
#
euso_l1.o: euso_l1.cpp
	gcc -c euso_l1.cpp `root-config --libs --cflags` -I `root-config --incdir`  -lstdc++
#
clean:
	rm euso_l1
	rm *.o
