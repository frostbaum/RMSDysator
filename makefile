FC=gfortran
FCFLAGS+= -O2
LDFLAGS+= -llapack -lblas
p_NAME := RMSDysator_prog

.PHONY: all, clean, distclean

all: $(p_NAME)

$(p_NAME): v3d_func_rep.o mk.o structure_c.o llist_c.o rmsd_m.o qsort_struc_m.o main.o
	$(FC) $(FCFLAGS) v3d_func_rep.o mk.o structure_c.o llist_c.o rmsd_m.o qsort_struc_m.o main.o $(LDFLAGS) -o $(p_NAME)

v3d_func_rep.o v3d_func_rep.mod: 3d-vectorOPs/v3d_func_rep.f95
	$(FC) $(FCFLAGS) -c 3d-vectorOPs/v3d_func_rep.f95
	@ touch v3d_func_rep.mod

mk.o mk.mod: MK-alg/mk.f95
	$(FC) $(FCFLAGS) -c MK-alg/mk.f95
	@ touch mk.mod

structure_c.o structure_c.mod: structure/structure_c.f95 v3d_func_rep.mod
	$(FC) $(FCFLAGS) -c structure/structure_c.f95
	@ touch structure_c.mod

llist_c.o llist_c.mod: llist_gen/llist_c.f95
	$(FC) $(FCFLAGS) -c llist_gen/llist_c.f95
	@ touch llist_c.mod

rmsd_m.o rmsd_m.mod: rmsd_m.f95 v3d_func_rep.mod mk.mod structure_c.mod llist_c.mod
	$(FC) $(FCFLAGS) -c rmsd_m.f95
	@ touch rmsd_m.mod

qsort_struc_m.o qsort_struc_m.mod: qsort_struc_m.f95 structure_c.mod
	$(FC) $(FCFLAGS) -c qsort_struc_m.f95
	@ touch qsort_struc_m.mod

main.o: main.f95 structure_c.mod qsort_struc_m.mod rmsd_m.mod
	$(FC) $(FCFLAGS) -c main.f95


#testprog.o: testprog.f95 rmsd_m.mod structure_c.mod
#	$(FC) $(FCFLAGS) -c testprog.f95

clean:
	@- $(RM) $(p_NAME)
	@- $(RM) *.o *.mod

distclean: clean
