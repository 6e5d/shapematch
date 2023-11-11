#include <stddef.h>

#include "../include/particle.h"
#include "../../modelobj/include/modelobj.h"

typedef struct {
	size_t plen;
	ShapematchParticle* ps;
} Shapematch;

void shapematch_cmass(Shapematch* sm, vec3 cmass);
void shapematch_init(Shapematch* sm, Modelobj* model);
void shapematch_step(Shapematch* sm);
void shapematch_deinit(Shapematch* sm);
