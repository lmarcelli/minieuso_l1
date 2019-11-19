#################################################################
####   Makefile to compile Level1 Mini-EUSO Software   ##########
#################################################################
euso_l1_tle: euso_l1_tle.o
	gcc -o euso_l1_tle euso_l1_tle.o `root-config --libs --cflags` -I `root-config --incdir`  -lstdc++
#
euso_l1_tle.o: euso_l1_tle.cpp
	gcc -c euso_l1_tle.cpp `root-config --libs --cflags` -I `root-config --incdir`  -lstdc++
#
clean:
	rm euso_l1_tle
	rm *.o
