all:
	g++ -std=c++11 -O3 -g -o meet_point meeting.cpp
	g++ -std=c++11 -O3 -o shortest shortest.cpp
clean:
	rm shortest meet_point
