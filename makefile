F90=mpiifort

files=globvars.f90 eigen.f90 auxfun.f90 symmetry.f90 matelem.f90  \
MPImat.f90 gvmopt.f90 svmopt.f90 io.f90 main.f90
objs=globvars.o eigen.o auxfun.o symmetry.o matelem.o  \
MPImat.o gvmopt.o svmopt.o IO.o main.o
mods=globvars.mod eigen.mod auxfun.mod symmetry.mod matelem.mod  \
mpimat.mod gvmopt.mod svmopt.mod io.mod 

LIBS1=-L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_scalapack_lp64
LIBS2=-L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_lapack95_lp64  
LIBS3=-L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_blas95_lp64
LIBS4=-L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_blacs_intelmpi_lp64
LIBS5=-L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_intel_lp64
LIBS6=-L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_sequential
LIBS7=-L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_core

main : $(objs)
	$(F90) -o main $(LIBS1) $(LIBS2) $(LIBS3) $(LIBS4) $(LIBS5) $(LIBS6) $(LIBS7) $(objs) 

main.o : main.f90
	$(F90) -c $?

globvars.o : globvars.f90
	$(F90) -c  $?

eigen.o :  eigen.f90
	$(F90) -c $? 

auxfun.o : auxfun.f90
	$(F90) -c $?

symmetry.o : symmetry.f90
	$(F90) -c $?

matelem.o : matelem.f90 
	$(F90) -c $?

MPImat.o : MPImat.f90
	$(F90) -c $?

svmopt.o : svmopt.f90 
	$(F90) -c $?

gvmopt.o : gvmopt.f90 
	$(F90) -c $?

IO.o : IO.f90
	$(F90) -c $?

clean :
	rm  $(objs) $(mods)	