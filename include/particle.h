#ifndef INCLUDEGUARD_SHAPEMATCH_PARTICLEH
#define INCLUDEGUARD_SHAPEMATCH_PARTICLEH

#include <cglm/cglm.h>

typedef struct {
	float* pos;
	vec3 pps;
	vec3 r0;
	float mass;
} ShapematchParticle;

#endif
