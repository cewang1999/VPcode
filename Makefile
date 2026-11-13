# 工作站 (Linux, g++) 编译: make
# 不用 OpenMP: make OMP=
CXX      ?= g++
OMP      ?= -fopenmp
CXXFLAGS ?= -O3 -std=c++17
# 写入输出文件头的代码版本 (无 git 时为 unknown)
GIT      := $(shell git describe --always --dirty 2>/dev/null || echo unknown)

vp: *.cpp *.h
	$(CXX) $(CXXFLAGS) $(OMP) -DVP_GIT='"$(GIT)"' *.cpp -o vp

clean:
	rm -f vp
