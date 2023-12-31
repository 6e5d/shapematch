#include <cglm/cglm.h>
#include <math.h>

#include "../../cglmh/build/cglmh.h"
#include "../../linalg/build/linalg.h"
#include "../../modelobj/include/modelobj.h"
#include "../include/shapematch.h"

const static size_t warn_length = 1000;

static void matsqrt(mat3 cov, mat3 sqr) {
	mat3 t, vec;
	vec3 val;
	linalg(eigen)((float*)cov, (float*)vec, (float*)val);
	if (!(val[0] >= 0.0f && val[1] >= 0.0f && val[2] >= 0.0f)) {
		printf("error: bad eigenvalue of positive semifinite\n");
		com_6e5d_cglmh_debug_vec3(val);
		return;
	}
	// vec * sqrt * t=inv
	glm_mat3_inv(vec, t);
	sqr[0][0] = sqrtf(val[0]);
	sqr[1][1] = sqrtf(val[1]);
	sqr[2][2] = sqrtf(val[2]);
	glm_mat3_mul(t, sqr, sqr);
	glm_mat3_mul(sqr, vec, sqr);
}

static void step1(Shapematch()* sm) {
	for (size_t idx = 0; idx < sm->plen; idx += 1) {
		Shapematch(Particle)* p = &sm->ps[idx];
		if (p->pos[1] < 0.0f) {
			p->pos[1] = 0.0f;
			p->pps[0] = p->pos[0];
			p->pps[2] = p->pos[2];
		}
		vec3 pps;
		glm_vec3_copy(p->pos, pps);
		vec3 dp;
		glm_vec3_sub(pps, p->pps, dp);
		dp[1] -= 0.0005f;
		glm_vec3_scale(dp, 0.999f, dp);
		glm_vec3_add(p->pos, dp, p->pos);
		glm_vec3_copy(pps, p->pps);
	}
}

void shapematch(step)(Shapematch()* sm) {
	step1(sm);
	vec3 cmass = {0};
	shapematch(cmass)(sm, cmass);
	mat3 cov = {0};
	for (size_t idx = 0; idx < sm->plen; idx += 1) {
		Shapematch(Particle)* p = &sm->ps[idx];
		vec3 r1 = {0};
		glm_vec3_sub(p->pos, cmass, r1);
		float a = r1[0]; float b = r1[1]; float c = r1[2];
		float x = p->r0[0]; float y = p->r0[1]; float z = p->r0[2];
		mat3 dcov = {
			{a * x, b * x, c * x},
			{a * y, b * y, c * y},
			{a * z, b * z, c * z},
		};
		com_6e5d_cglmh_mat3_add(cov, dcov, cov);
	}
	mat3 t, sqr = {0};
	glm_mat3_transpose_to(cov, t);
	glm_mat3_mul(t, cov, t);
	// sqr is the symmetric part to be excluded
	matsqrt(t, sqr);
	glm_mat3_inv(sqr, t);
	glm_mat3_mul(cov, t, t);
	for (size_t idx = 0; idx < sm->plen; idx += 1) {
		Shapematch(Particle)* p = &sm->ps[idx];
		vec3 dp;
		glm_mat3_mulv(t, p->r0, dp);
		glm_vec3_add(dp, cmass, dp);
		glm_vec3_sub(dp, p->pos, dp);
		// rigidity compliance
		glm_vec3_scale(dp, 0.6f, dp);
		glm_vec3_add(p->pos, dp, p->pos);
	}
}

void shapematch(cmass)(Shapematch()* sm, vec3 cmass) {
	glm_vec3_zero(cmass);
	for (size_t idx = 0; idx < sm->plen; idx += 1) {
		Shapematch(Particle)* p = &sm->ps[idx];
		glm_vec3_add(p->pos, cmass, cmass);
	}
	glm_vec3_scale(cmass, 1.0f / (float)sm->plen, cmass);
}

void shapematch(init)(Shapematch()* sm, Modelobj()* model) {
	if (model->v_len > warn_length) {
		// though cmass/cov precision can be solved kahan addition,
		// for physical model we generally avoid to have many points.
		// instead of solving precision, user should
		// consider separate render and physical model
		// and determine the rest vertex positions by interpolation
		printf("model contains too many vertices %zu>%zu\n",
			model->v_len, warn_length);
	}
	sm->plen = model->v_len;
	sm->ps = calloc(model->v_len, sizeof(Shapematch(Particle)));
	for (size_t idx = 0; idx < model->v_len; idx += 1) {
		Shapematch(Particle)* p = &sm->ps[idx];
		p->mass = 1.0f;
		p->pos = model->vs[idx];
		glm_vec3_copy(p->pos, p->pps);
	}
	vec3 cmass = {0};
	shapematch(cmass)(sm, cmass);
	for (size_t idx = 0; idx < model->v_len; idx += 1) {
		Shapematch(Particle)* p = &sm->ps[idx];
		glm_vec3_sub(p->pos, cmass, p->r0);
	}
}

void shapematch(deinit)(Shapematch()* sm) {
	free(sm->ps);
}
