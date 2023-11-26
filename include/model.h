#ifndef INCLUDEGUARD_SHAPEMATCH_MODELH
#define INCLUDEGUARD_SHAPEMATCH_MODELH

#include <stddef.h>

#include "../../modelobj/include/modelobj.h"
#include "../include/particle.h"

typedef struct {
size_t plen;
ShapematchParticle* ps;
} Shapematch;

void shapematch_cmass(Shapematch* sm, vec3 cmass);
void shapematch_init(Shapematch* sm, Modelobj* model);
void shapematch_step(Shapematch* sm);
void shapematch_deinit(Shapematch* sm);

#endif
