all: desert_finder desert_comparison filter_bed_length grep_bed_region

desert_finder.o: desert_finder.cpp utils.cpp
	g++ -std=c++11 -Wall -c desert_finder.cpp -o desert_finder.o

desert_comparison.o: desert_comparison.cpp utils.cpp
	g++ -std=c++11 -Wall -c desert_comparison.cpp -o desert_comparison.o

filter_bed_length.o: filter_bed_length.cpp
	g++ -std=c++11 -Wall -c filter_bed_length.cpp -o filter_bed_length.o

grep_bed_region.o: grep_bed_region.cpp
	g++ -std=c++11 -Wall -c grep_bed_region.cpp -o grep_bed_region.o

utils.o: utils.cpp
	g++ -std=c++11 -Wall -c utils.cpp -o utils.o

desert_finder: desert_finder.o utils.o
	g++ desert_finder.o utils.o -o desert_finder

desert_comparison: desert_comparison.o utils.o
	g++ desert_comparison.o utils.o -o desert_comparison

filter_bed_length: filter_bed_length.o
	g++ filter_bed_length.o -o filter_bed_length

grep_bed_region: grep_bed_region.o
	g++ grep_bed_region.o -o grep_bed_region

clean:
	rm *.o desert_finder desert_comparison filter_bed_length grep_bed_region
