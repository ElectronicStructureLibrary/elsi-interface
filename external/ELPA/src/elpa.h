#ifndef ELPA_H
#define ELPA_H

#include <limits.h>
#include <complex.h>

#include "elpa_version.h"

struct elpa_struct;
typedef struct elpa_struct *elpa_t;

struct elpa_autotune_struct;
typedef struct elpa_autotune_struct *elpa_autotune_t;

#include "elpa_constants.h"
#include "elpa_generated.h"
#include "elpa_generic.h"

const char *elpa_strerr(int elpa_error);

#endif
