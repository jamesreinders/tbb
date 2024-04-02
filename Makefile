all: pi
	./pi

pi: pi.cpp
	icpx -tbb -o pi pi.cpp
