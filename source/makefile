F90=mpiifort

objs=globvars.o eigen.o auxfun.o symmetry.o matelem.o  \
MPImat.o gvmopt.o svmopt.o IO.o main.o

libs=-L/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_scalapack_lp64 \
-lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_blacs_intelmpi_lp64 \
-lmkl_intel_lp64 -lmkl_sequential -lmkl_core

main : $(objs)
	$(F90) -o main $(libs) $(objs) 

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
	rm -f *.o
	rm -f *.mod
