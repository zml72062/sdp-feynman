SRCDIR    = src
INCDIR    = include
OBJDIR    = build
CXX       = g++
CPPFLAGS  = -I./${INCDIR}
CXXFLAGS  = -O2
LDFLAGS   = -lyaml-cpp -lcln -lginac -lgsl -lgslcblas -lsdpa -ldmumps_seq -llapack -lblas

OBJS      = ${OBJDIR}/config.o \
			${OBJDIR}/utils.o \
			${OBJDIR}/ibp.o \
			${OBJDIR}/parse.o \
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
	rm -f master

