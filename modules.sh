module load openmpi/4.1.1
export PATH=$HOME/julia-1.6.3/bin:$PATH
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/apps/openmpi/4.1.1
export JULIA_PETSC_LIBRARY="$HOME/petsc-install/lib/libpetsc"
export P4EST_ROOT_DIR=$HOME/p4est-install/
export UCX_ERROR_SIGNALS="SIGILL,SIGBUS,SIGFPE"
export GRIDAP_PARDISO_LIBGOMP_DIR="/half-root/lib/gcc/x86_64-redhat-linux/8/"
