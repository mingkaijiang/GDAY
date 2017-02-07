#ifndef WRITE_OUT_H
#define WRITE_OUT_H


#include "gday.h"
#include "utilities.h"

void  open_output_file(control *, char *, FILE **);
void  write_output_header(control *, params *, FILE **);
void  write_annual_outputs_ascii(control *, fluxes *, state *, int, int);
int   write_final_state(control *, params *p, state *);
int   ohandler(char *, char *, char *, control *, params *p, state *, int *);


#endif /* WRITE_OUT_H */
