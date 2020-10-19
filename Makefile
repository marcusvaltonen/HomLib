.PHONY: test clean example lint

CXX=g++
RM=rm -f
CPPFLAGS=--std=c++11 -O2
CPPTESTFLAGS=--coverage -g -O0
LDTESTFLAGS=--coverage
INCLUDES=\
	-I./includes/HomLib \
	-I./src/helpers \
	-I./src/solvers/fitzgibbon_cvpr_2001 \
	-I./src/solvers/kukelova_cvpr_2015 \
	-I./src/solvers/valtonenornhag_arxiv_2020a \
	-I./src/solvers/valtonenornhag_arxiv_2020b

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	INCLUDES += -I/usr/include/eigen3
endif
ifeq ($(UNAME_S),Darwin)
	INCLUDES += -I/usr/local/Cellar/eigen/3.3.7/include
endif

SRCS=\
	src/helpers/generate_problem_instance.cpp \
	src/helpers/gj.cpp \
	src/helpers/normalize2dpts.cpp \
	src/helpers/radial.cpp \
	src/helpers/roots.cpp \
	src/solvers/fitzgibbon_cvpr_2001/get_fitzgibbon_cvpr_2001.cpp \
	src/solvers/kukelova_cvpr_2015/get_kukelova_cvpr_2015.cpp \
	src/solvers/kukelova_cvpr_2015/solver_kukelova_cvpr_2015.cpp \
	src/solvers/valtonenornhag_arxiv_2020a/get_valtonenornhag_arxiv_2020a_fHf.cpp \
	src/solvers/valtonenornhag_arxiv_2020a/solver_valtonenornhag_arxiv_2020a_fHf.cpp \
	src/solvers/valtonenornhag_arxiv_2020b/get_valtonenornhag_arxiv_2020b_fHf.cpp \
	src/solvers/valtonenornhag_arxiv_2020b/get_valtonenornhag_arxiv_2020b_frHfr.cpp \
	src/solvers/valtonenornhag_arxiv_2020b/solver_valtonenornhag_arxiv_2020b_fHf.cpp \
	src/solvers/valtonenornhag_arxiv_2020b/solver_valtonenornhag_arxiv_2020b_frHfr.cpp

EXAMPLE_SRCS=example.cpp

TEST_SRCS=\
	tests/test_fitzgibbon_cvpr_2001.cpp \
	tests/test_helpers.cpp \
	tests/test_kukelova_cvpr_2015.cpp \
	tests/test_main.cpp \
	tests/test_valtonenornhag_arxiv_2020a.cpp \
	tests/test_valtonenornhag_arxiv_2020b.cpp

TEST_INCLUDES=\
	-I./tests

OBJS = $(patsubst %.cpp,%.o,$(SRCS))
EXAMPLE_OBJS = $(patsubst %.cpp,%.o,$(EXAMPLE_SRCS))
TEST_OBJS = $(patsubst %.cpp,%.o,$(TEST_SRCS))

all: example

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(TEST_INCLUDES) -c -o $@ $<

test: CPPFLAGS += $(CPPTESTFLAGS)
test: CPPFLAGS:=$(filter-out -O2,$(CPPFLAGS))
test: $(OBJS) $(TEST_OBJS)
	$(CXX) -o run_tests $(OBJS) $(TEST_OBJS) $(LDTESTFLAGS)
	@./run_tests

example: $(OBJS) $(EXAMPLE_OBJS)
	$(CXX) -o example $(OBJS) $(EXAMPLE_OBJS)

lint:
	cpplint --linelength=119 --verbose=1 --exclude=lib --exclude=tests/catch2 --recursive --quiet .

clean:
	$(RM) $(OBJS)
	$(RM) $(TEST_OBJS)
	$(RM) $(EXAMPLE_OBJS)
	$(RM) $(patsubst %.cpp,%.gcno,$(EXAMPLE_SRCS))
	$(RM) $(patsubst %.cpp,%.gcno,$(TEST_SRCS))
	$(RM) $(patsubst %.cpp,%.gcda,$(TEST_SRCS))
	$(RM) $(patsubst %.cpp,%.gcno,$(SRCS))
	$(RM) $(patsubst %.cpp,%.gcda,$(SRCS))
	$(RM) example
	$(RM) run_tests
