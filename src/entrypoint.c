#include <R.h>
#include <Rinternals.h>

void R_init_ptw_rust_extendr(DllInfo *dll);

void R_init_ptwRust(DllInfo *dll) {
    R_init_ptw_rust_extendr(dll);
}

