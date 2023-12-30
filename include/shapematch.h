#include <cglm/cglm.h>

#include "../../modelobj/include/modelobj.h"

typedef struct {
	float* pos;
	vec3 pps;
	vec3 r0;
	float mass;
} Shapematch(Particle);

typedef struct {
	size_t plen;
	Shapematch(Particle)* ps;
} Shapematch();

void shapematch(cmass)(Shapematch()* sm, vec3 cmass);
void shapematch(init)(Shapematch()* sm, Modelobj()* model);
void shapematch(step)(Shapematch()* sm);
void shapematch(deinit)(Shapematch()* sm);
