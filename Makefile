NO_SDPA_LIB = false

SRCDIR    = src
INCDIR    = include
OBJDIR    = build
CXX       = g++

CPPFLAGS  = -I./${INCDIR}
ifeq (${NO_SDPA_LIB}, true)
CPPFLAGS  += -DNO_SDPA_LIB
endif
ifeq (${NO_GSL}, true)
CPPFLAGS  += -DNO_GSL
endif

CXXFLAGS  = -O2

LDFLAGS   = -lyaml-cpp -lcln -lginac
ifneq (${NO_SDPA_LIB}, true)
LDFLAGS   += -lsdpa -ldmumps_seq -llapack -lblas
endif
ifneq (${NO_GSL}, true)
LDFLAGS   += -lgsl -lgslcblas
endif

OBJS      = ${OBJDIR}/config.o \
			${OBJDIR}/utils.o \
			${OBJDIR}/ibp.o \
			${OBJDIR}/cache.o \
			${OBJDIR}/subprocess.o \
			${OBJDIR}/parse.o \
			${OBJDIR}/dimshift.o \
			${OBJDIR}/generate.o \
			${OBJDIR}/sdpa.o \
			${OBJDIR}/solver.o \
			${OBJDIR}/main.o

all: pre master

master: ${OBJS}
	${CXX} ${OBJS} -o master ${LDFLAGS}

${OBJS}: ${OBJDIR}/%.o: ${SRCDIR}/%.cpp
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $^ -o $@

.PHONY: pre
pre:
	mkdir -p ${OBJDIR}

.PHONY: clean
clean:
	rm -rf ${OBJDIR}
	rm -rf tmp
	rm -f master

