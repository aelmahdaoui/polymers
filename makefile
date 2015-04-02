FC = gfortran
FFLAGS =  -march=native -O3 -Wall -Wextra -Wtabs -fcheck=all  -mno-avx -ffast-math -fopenmp

LDFLAGS = -fopenmp
LIBS = -llapack -lblas

FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)




OBJS = 
OBJS += polymers_init.o
OBJS += polymer_dynamics.o
OBJS += polymers_results.o
OBJS += polymers.o



all: polymers



polymers: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<


.PHONY: clean
clean:
	$(RM) polymers  $(OBJS) *.mod 
	rm -rf obj/ mod/





