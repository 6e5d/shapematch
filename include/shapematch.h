#ifndef INCLUDEGUARD_SHAPEMATCH
#define INCLUDEGUARD_SHAPEMATCH

#include <cglm/cglm.h>

#include "../../modelobj/include/modelobj.h"
#define CglmMat3 mat3
#define CglmVec3 vec3

typedef struct {
	float* pos;
	CglmVec3 pps;
	CglmVec3 r0;
	float mass;
} Shapematch(Particle);

typedef struct {
	size_t plen;
	Shapematch(Particle)* ps;
} Shapematch();

void shapematch(cmass)(Shapematch()* sm, CglmVec3 cmass);
void shapematch(init)(Shapematch()* sm, Modelobj()* model);
void shapematch(step)(Shapematch()* sm);
void shapematch(deinit)(Shapematch()* sm);

#endif
