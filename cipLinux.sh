mpiifort "$1".f90 -o "$1" -O3  -I"${MKLROOT}/include" -L/usr/local/lib/libcfitsio.a ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -lcfitsio -Wl,--end-group -liomp5 -lpthread -lm -ldl

