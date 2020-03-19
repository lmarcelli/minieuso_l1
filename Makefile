#################################################################
####   Makefile to compile Level1 Mini-EUSO Software   ##########
#################################################################
euso_l1_tle_v7: euso_l1_tle_v7.o
	gcc -o euso_l1_tle_v7 euso_l1_tle_v7.o `root-config --libs --cflags` -I `root-config --incdir`  -lstdc++
#
euso_l1_tle_v7.o: euso_l1_tle_v7.cpp
	gcc -c euso_l1_tle_v7.cpp `root-config --libs --cflags` -I `root-config --incdir`  -lstdc++
#
clean:
	rm euso_l1_tle_v7
	rm *.o
