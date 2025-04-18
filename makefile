PRES_DIR = .
WORK_DIR = ./testing
console_DIR = ./console
# DIRR = ${PRES_DIR}/randomc
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

# 添加头文件搜索路径
HEADERS := $(wildcard ${DIRG}/*.h) $(wildcard ${DIRS}/*.h)

# 添加并行编译支持
MAKEFLAGS += -j8

default: manpower

srun: ${EXE}
	@echo '______________________________________________ sbatch ______________________________________________'
	@echo ''
	@qsub ${console_DIR}/slurm.bscc.t6.txt
	@echo '____________________________________________________________________________________________________'
	@echo ''

qsub: ${norg}
	@echo '_______________________________________________ qsub _______________________________________________'
	@echo ''
	@qsub ${WORK_DIR}/qsub.DPC++CPU.mpi.txt > ${WORK_DIR}/jmwang.job.txt
	@chmod +x ${WORK_DIR}/job.process.DPC++CPU.txt
	@${WORK_DIR}/job.process.DPC++CPU.txt
	@rm -rf ${WORK_DIR}/jmwang.job.txt
	@echo '____________________________________________________________________________________________________'
	@echo ''

${norg}: compile $(DEPS) $(OBJS)
	@echo '_______________________________________________ link _______________________________________________'
	@echo ''
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -o ${norg} $(OBJS) $(LIBS)

manpower: ${norg}
	@echo '_____________________________ Now you can test by: mpirun -n 24 ./inorg _____________________________'

start:
	@echo '______________________________________________ start ______________________________________________'
	@echo ''

explain:
	@echo '_____________________________________________ explain _____________________________________________'
	@echo ''
	@echo "the following information represents your prgram"
	@echo "final norg name: $(norg)"
	@echo "source files: $(SRCS)"
	@echo "object files: $(OBJS)"
	@echo "dependency files: $(DEPS)"

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
	# -rm -rf bi/*.txt
	# -rm -rf bi/*.out
	# -rm -rf io/output.*
	# -rm -rf io/*.txt
	# -rm -rf tso/*
	# -rm -rf jmwang.*
	# -rm -rf *jmw*
	# -rm -rf log.*
deepclean:
	# -rm -rf general/*.cpp
	# -rm -rf general/*.h
	# -rm -rf special/*.cpp
	# -rm -rf special/*.h
	# -rm -rf randomc/*.cpp
	# -rm -rf randomc/*.h
	

depend:$(DEPS)
	@echo "dependencies are now up-to-date."

-include $(DEPS)
