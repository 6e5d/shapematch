#ifndef INCLUDEGUARD_SHAPEMATCH
#define INCLUDEGUARD_SHAPEMATCH

#include <cglm/cglm.h>

#include "../../modelobj/include/modelobj.h"

typedef struct {
	float* pos;
	vec3 pps;
	vec3 r0;
	float mass;
} ShapematchParticle;

typedef struct {
	size_t plen;
	ShapematchParticle* ps;
} Shapematch;

void shapematch_cmass(Shapematch* sm, vec3 cmass);
void shapematch_init(Shapematch* sm, Modelobj* model);
void shapematch_step(Shapematch* sm);
void shapematch_deinit(Shapematch* sm);

#endif
