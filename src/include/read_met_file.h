#ifndef READ_MET_H
#define READ_MET_H

#include <stdio.h>
#include <stdlib.h>

#include "gday.h"
#include "utilities.h"
void    read_daily_met_data_simple(char **, control *, met_arrays *, params *);

void    read_daily_met_data(char **, control *, met_arrays *);

#endif /* READ_MET_H */
