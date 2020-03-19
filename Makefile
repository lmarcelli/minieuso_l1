#################################################################
####   Makefile to compile Level1 Mini-EUSO Software   ##########
#################################################################
euso_l1_tle_pu_v7r1: euso_l1_tle_pu_v7r1.o
	gcc -o euso_l1_tle_pu_v7r1 euso_l1_tle_pu_v7r1.o `root-config --libs --cflags` -I `root-config --incdir`  -lstdc++
#
euso_l1_tle_pu_v7r1.o: euso_l1_tle_pu_v7r1.cpp
	gcc -c euso_l1_tle_pu_v7r1.cpp `root-config --libs --cflags` -I `root-config --incdir`  -lstdc++
#
clean:
	rm euso_l1_tle_pu_v7r1
	rm *.o
