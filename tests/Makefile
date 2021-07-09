CXXFLAGS = -std=c++11 -Wall -Wunused-result -O3 -DNDEBUG -g -ffast-math

# !! Error w/ valgrind -- does it actually improve performance?
# CXXFLAGS += -std=c++11 -mtune=native -march=native -Wall -Wunused-result -O3 -DNDEBUG -g -fopenmp -ffast-math

APP=test

all: $(APP)

$(APP): $(APP).cpp
	g++ $(CXXFLAGS) -o $(APP) $(APP).cpp ../lap/lapmod.cpp -I../lap

clean:
	rm -f $(APP)
