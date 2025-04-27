PRES_DIR = .
WORK_DIR = .
DIRG = ${PRES_DIR}/src/gen
DIRS = ${PRES_DIR}/src
norg  = ${WORK_DIR}/inorg
VPATH = ${DIRG}:${DIRS}

CC       = mpiicpc
# CPPFLAGS = -std=c++17 -qmkl -O3
CPPFLAGS = -std=c++17 -mkl -O3
CFLAGS   =
CXXFLAGS = $(CFLAGS)
COMPILE  = $(CC) $(CPPFLAGS) $(CXXFLAGS) -c

SRCS := $(wildcard ${DIRG}/*.cpp) $(wildcard ${DIRS}/*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))
DEPS := $(patsubst %.cpp,%.d,$(SRCS))

# Add header file search paths
HEADERS := $(wildcard ${DIRG}/*.h) $(wildcard ${DIRS}/*.h)

# Add parallel compilation support
MAKEFLAGS += -j8

default: manpower


${norg}: compile $(DEPS) $(OBJS)
	@echo '_______________________________________________ link _______________________________________________'
	@echo ''
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -o ${norg} $(OBJS) $(LIBS)

manpower: ${norg}
	@echo '_____________________________ Now you can test by: mpirun -n 6 ./inorg _____________________________'

compile:
	@echo '_____________________________________________ compile _____________________________________________'
	@echo ''

# Optimize dependency file generation rules
%.d: %.cpp
	@echo "Generating dependency for $<"
	@$(CC) -MM -MP -MF"$@" -MT"$*.o $@" $(CPPFLAGS) $< > $@
	@if [ -f $*.h ]; then \
		echo "$*.o: $*.h" >> $@; \
		echo "$@: $*.h" >> $@; \
	fi

# Optimize object file compilation rules
%.o: %.cpp
	@echo "Compiling $<"
	@$(COMPILE) -o $@ $<

clean:
	-rm -rf $(OBJS) $(DEPS) $(norg)

clear:
	-rm -rf testing/edmft_back_up
	-rm -rf testing/g0imp.txt
	-rm -rf testing/Gf.out
	-rm -rf testing/h0.txt
	-rm -rf testing/hb_fit.txt
	-rm -rf testing/hb_read.txt
	-rm -rf testing/inorg
	-rm -rf testing/log.norg
	-rm -rf testing/nmat.txt
	-rm -rf testing/nohup.txt
	-rm -rf testing/ose_hop
	-rm -rf testing/Sig.out
	

depend:$(DEPS)
	@echo "dependencies are now up-to-date."

-include $(DEPS)
