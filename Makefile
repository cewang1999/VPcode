# 工作站 (Linux, g++) 编译: make
# 不用 OpenMP: make OMP=
CXX      ?= g++
OMP      ?= -fopenmp
CXXFLAGS ?= -O3 -std=c++17

vp: *.cpp *.h
	$(CXX) $(CXXFLAGS) $(OMP) *.cpp -o vp

clean:
	rm -f vp
