/*
 * evolsim.c
 *
 *  Created on: Jul 26, 2016
 *      Author: root
 */

#include "evolsim.h"
#define GLEW_STATIC
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include <GLFW/glfw3.h>
#include "font.h"
#include <errno.h>
#include "xstring.h"
#include "globals.h"
#include <stdio.h>
#include "gui.h"
#include <math.h>
#include "smem.h"
#include "pqueue.h"
#include "wincompat.h"

void rect(float x, float y, float width, float height) {
	glBegin (GL_QUADS);
	glVertex2f(x, y);
	glVertex2f(x, y + height);
	glVertex2f(x + width, y + height);
	glVertex2f(x + width, y);
	glEnd();
}

void circle(float x, float y, float rad) {
	glBegin (GL_TRIANGLE_FAN);
	glVertex2f(x, y);
	for (int i = 20; i >= 0; i--) {
		glVertex2f(x + (rad * cos((float) i * 2. * M_PI / 20.)), y + (rad * sin((float) i * 2. * M_PI / 20.)));
	}
	glEnd();
}

void line(float x1, float y1, float x2, float y2) {
	glBegin (GL_LINES);
	glVertex2f(x1, y1);
	glVertex2f(x2, y2);
	glEnd();
}

void text_center(char* str, int x, int y, float scale) {
	glEnable (GL_TEXTURE_2D);
	glEnable (GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPushMatrix();
	glTranslatef(x, y, 0.);
	glScalef(scale, scale, scale);
	drawString(str, -((float) stringWidth(str) / 2.), -10., 0);
	glPopMatrix();
	glDisable(GL_BLEND);
	glDisable(GL_TEXTURE_2D);
}

void text_right(char* str, int x, int y, float scale) {
	glEnable (GL_TEXTURE_2D);
	glEnable (GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPushMatrix();
	glTranslatef(x, y, 0.);
	glScalef(scale, scale, scale);
	drawString(str, -((float) stringWidth(str)), -10., 0);
	glPopMatrix();
	glDisable(GL_BLEND);
	glDisable(GL_TEXTURE_2D);
}

void text_left(char* str, int x, int y, float scale) {
	glEnable (GL_TEXTURE_2D);
	glEnable (GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPushMatrix();
	glTranslatef(x, y, 0.);
	glScalef(scale, scale, scale);
	drawString(str, 0., -10., 0);
	glPopMatrix();
	glDisable(GL_BLEND);
	glDisable(GL_TEXTURE_2D);
}

struct creature* cpop[1000];
struct creature* ccsimc[1000];

void evol_render(float partialTick) {
	glViewport(0, 0, width, height);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0., width, height, 0., 1000., 3000.);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0., 0., -2000.);
	glClearColor(1., 1., 1., 1.);
	glColor4f(100. / 255., 200. / 255., 100 / 255., 1.);
	rect(width / 2. - 200., 300., 400., 200.);
	glColor4f(1., 1., 1., 1.);
	text_center("EVOLUTION!", width / 2, 200, 6.);
	text_center("START", width / 2, 430, 6.);
}

void evol_button(int action, int mods, double x, double y) {
	guistate_set(1);
}

void evol_init() {
	if (loadFont("assets" PATH_SEP "default_font.png", "assets" PATH_SEP "default_fontWidth.bin", TX_FONT)) {
		printf("Error loading font: %s\n", strerror(errno));
	}
	setFont (TX_FONT);
	struct gui_button btn;
	btn.action = 0;
	btn.width = (int) (400.);
	btn.height = (int) (200.);
	btn.x = (int) (width / 2. - 200.);
	btn.y = (int) (300.);
	gui_addbutton(0, btn);
}

int mkcr_gensel = -1;
int mkcr_gen = 0;

struct histogram {
		float* data;
		size_t data_count;
		float median;
};

void drawHistogram(struct histogram* hg, int x, int y, int hw, int hh) {
	int maxH = 1;
	for (size_t i = 0; i < hg->data_count; i++) {
		if (hg->data[i] > maxH) {
			maxH = hg->data[i];
		}
	}
	glColor4f(200. / 255., 200. / 255., 200. / 255., 1.);
	rect(x, y, hw, hh);
	float barW = (float) hw / barLen;
	float multiplier = (float) hh / maxH * 0.9;
	//stroke(128);
	//strokeWeight(2);
	//glLineWidth(128.);
	int unit = 100;
	if (maxH < 300) unit = 50;
	if (maxH < 100) unit = 20;
	if (maxH < 50) unit = 10;
	char pbuf[256];
	for (int i = 0; i < hh / multiplier; i += unit) {
		float theY = y + hh - i * multiplier;
		line(x, theY, x + hw, theY);
		if (i == 0) theY -= 5;
		snprintf(pbuf, 256, "%i", i);
		text_left(pbuf, x + hw + 5, theY + 7, 2.);
	}
	for (int i = minBar; i <= maxBar; i += 10) {
		/*if (i == 0) {
		 stroke(0, 0, 255);
		 } else {
		 stroke(128);
		 }*/
		float theX = x + (i - minBar) * barW;
		snprintf(pbuf, 256, "%.1f", (float) i / histBarsPerMeter);
		text_center(pbuf, theX, y + hh + 14, 2.);
		line(theX, y, theX, y + hh);
	}
	for (int i = 0; i < hg->data_count; i++) {
		float h = hg->data[i] * multiplier;
		if (hh < h) h = hh;
		if (i + minBar == (int) hg->median) {
			glColor4f(1., 0., 0., 1.);
		} else {
			glColor4f(0., 0., 0., 1.);
		}
		rect(x + i * barW, y + hh - h, barW, h);
	}
}

struct graph {
		float** data;
		size_t data_count;
		size_t data_subcount;
};

void drawGraphImage(struct graph* gr, int x, int y, int width, int height, float best, float worst) { // size 975, 570 @ 50, 180, downsized to 650, 380
	glColor4f(220. / 255., 220. / 255., 220. / 255., 1.);
	rect(x, y, width, height);
	if (gr->data_subcount >= 0) {
		float mh = (float) height / (best - worst);
		float gw = (float) width / (float) (gr->data_subcount == 0 ? 1 : gr->data_subcount);
		float zero = (best / (best - worst)) * 380.;
		float unit = 3. * log(best - worst) / log(10) - 2;
		if (fmod(unit + 90., 3.) < 1.) {
			unit = pow(10, (int) (unit / 3));
		} else if (fmod(unit + 90., 3.) < 2.) {
			unit = pow(10, (int) ((unit - 1) / 3)) * 2;
		} else {
			unit = pow(10, (int) ((unit - 2) / 3)) * 5;
		}
		glColor4f(150. / 255., 150. / 255., 150. / 255., 1.);
		char pbuf[256];
		for (float i = (int) (((worst - (best - worst) / 18.0) / unit) + .9999) * unit; i < best + (best - worst) / 18.0; i += unit) {
			float lineY = y - i * mh + zero;
			line(x, lineY, width + x, lineY);
			snprintf(pbuf, 256, "%.01f %s", i, fitnessUnit);
			text_right(pbuf, x - 5, lineY + 4, 2.);
		}
		for (int i = 0; i < gr->data_count; i++) { // 29
			int k;
			if (i == 28) {
				k = 14;
			} else if (i < 14) {
				k = i;
			} else {
				k = i + 1;
			}
			/*if (k == 14) {
			 graphImage.stroke(255, 0, 0, 255);
			 graphImage.strokeWeight(5);
			 } else {
			 stroke(0);
			 if (k == 0 || k == 28 || (k >= 10 && k <= 18)) {
			 graphImage.strokeWeight(3);
			 } else {
			 graphImage.strokeWeight(1);
			 }
			 }*/
			for (int j = 0; j < gr->data_subcount; j++) {
				line(x + j * gw, (-gr->data[k][j]) * mh + zero + y, x + (j + 1) * gw, (-gr->data[k][j + 1]) * mh + zero + y);
			}
		}
	}
}

struct histogram mkcr_mhg;
struct graph mkcr_mg;

void drawSegBars(int x, int y, int width, int height) {
	glColor4f(128. / 255., 128. / 255., 128. / 255., 1.);
	rect(x, y, width, height);
}

float evol_r() {
	return pow((drand48() * 2.) - 1., 19);
}

int evol_rInt() {
	return (int) ((drand48() * 2.) - 1.);
}

struct node {
		float x;
		float y;
		float vx;
		float vy;
		float px;
		float py;
		float pvx;
		float pvy;
		float m;
		float f;
		float value;
		float valueToBe;
		int op;
		int axon1;
		int axon2;
		int safeInput;
		float pressure;
};

struct muscle {
		int axon;
		int c1;
		int c2;
		float len;
		float rigidity;
		float previousTarget;
};

struct creature {
		struct node** nodes;
		size_t node_count;
		struct muscle** muscles;
		size_t muscle_count;
		float d;
		int id;
		int alive;
		float creatureTimer;
		float mutability;
};

struct sim_rect {
		int x1;
		int y1;
		int x2;
		int y2;
};

struct simulation {
		float totalNodeNausea;
		float averageNodeNausea;
		struct sim_rect* rects;
		size_t rect_count;
		size_t timer;
		float energy;
};

struct node* node_new(float tx, float ty, float tvx, float tvy, float tm, float tf, float val, int op, int axon1, int axon2) {
	struct node* n = smalloc(sizeof(struct node));
	n->px = n->x = tx;
	n->py = n->y = ty;
	n->pvx = n->vx = tvx;
	n->pvy = n->vy = tvy;
	n->m = tm;
	n->f = tf;
	n->value = n->valueToBe = val;
	n->op = op;
	n->axon1 = axon1;
	n->axon2 = axon2;
	n->pressure = 0;
	n->safeInput = 0;
	return n;
}

void node_applyForce(struct node* n, struct simulation* sim) {
	n->vx *= airFriction;
	n->vy *= airFriction;
	n->y += n->vy;
	n->x += n->vx;
	float dmx = n->vx - n->pvx;
	float dmy = n->vy - n->pvy;
	float accel = sqrt(dmx * dmx + dmy * dmy);
	if (sim != NULL) sim->totalNodeNausea += accel * accel * nauseaUnit;
	n->pvx = n->vx;
	n->pvy = n->vy;
}

void node_applyGravity(struct node* n) {
	n->vy += gravity;
}

void node_pressAgainstGround(struct node* n, float groundY) {
	float dif = n->y - (groundY - n->m / 2);
	n->pressure += dif * pressureUnit;
	n->y = (groundY - n->m / 2);
	n->vy = 0;
	n->x -= n->vx * n->f;
	if (n->vx > 0) {
		n->vx -= n->f * dif * friction;
		if (n->vx < 0) {
			n->vx = 0;
		}
	} else {
		n->vx += n->f * dif * friction;
		if (n->vx > 0) {
			n->vx = 0;
		}
	}
}

void node_hitWalls(struct node* n, struct simulation* sim) {
	n->pressure = 0.;
	float dif = n->y + n->m / 2;
	if (dif >= 0 && haveGround) {
		node_pressAgainstGround(n, 0);
	}
	if (n->y > n->py && hazelStairs >= 0) {
		float bottomPointNow = n->y + n->m / 2;
		float bottomPointPrev = n->py + n->m / 2;
		int levelNow = (int) (ceil(bottomPointNow / hazelStairs));
		int levelPrev = (int) (ceil(bottomPointPrev / hazelStairs));
		if (levelNow > levelPrev) {
			float groundLevel = levelPrev * hazelStairs;
			node_pressAgainstGround(n, groundLevel);
		}
	}
	for (size_t i = 0; i < sim->rect_count; i++) {
		struct sim_rect* r = &sim->rects[i];
		int flip = 0;
		float px = 0.;
		float py = 0.;
		if (abs(n->x - (r->x1 + r->x2) / 2) <= (r->x2 - r->x1 + n->m) / 2 && abs(n->y - (r->y1 + r->y2) / 2) <= (r->y2 - r->y1 + n->m) / 2) {
			if (n->x >= r->x1 && n->x < r->x2 && n->y >= r->y1 && n->y < r->y2) {
				float d1 = n->x - r->x1;
				float d2 = r->x2 - n->x;
				float d3 = n->y - r->y1;
				float d4 = r->y2 - n->y;
				if (d1 < d2 && d1 < d3 && d1 < d4) {
					px = r->x1;
					py = n->y;
				} else if (d2 < d3 && d2 < d4) {
					px = r->x2;
					py = n->y;
				} else if (d3 < d4) {
					px = n->x;
					py = r->y1;
				} else {
					px = n->x;
					py = r->y2;
				}
				flip = 1;
			} else {
				if (n->x < r->x1) {
					px = r->x1;
				} else if (n->x < r->x2) {
					px = n->x;
				} else {
					px = r->x2;
				}
				if (n->y < r->y1) {
					py = r->y1;
				} else if (n->y < r->y2) {
					py = n->y;
				} else {
					py = r->y2;
				}
			}
			float dx = n->x - px;
			float dy = n->y - py;
			float distance = sqrt(dx * dx + dy * dy);
			float rad = n->m / 2;
			float wallAngle = atan2(py - n->y, px - n->x);
			if (flip) {
				wallAngle += M_PI;
			}
			if (distance < rad || flip) {
				dif = rad - distance;
				n->pressure += dif * pressureUnit;
				float multi = rad / distance;
				if (flip) {
					multi = -multi;
				}
				n->x = (n->x - px) * multi + px;
				n->y = (n->y - py) * multi + py;
				float veloAngle = atan2(n->vy, n->vx);
				float veloMag = sqrt(n->vx * n->vx + n->vy * n->vy);
				float relAngle = veloAngle - wallAngle;
				float relY = sin(relAngle) * veloMag * dif * friction;
				n->vx = -sin(relAngle) * relY;
				n->vy = cos(relAngle) * relY;
			}
		}
	}
	n->py = n->y;
	n->px = n->x;
}

void node_doMath(struct node* n, struct simulation* sim, struct creature* cr) {
	float av1 = n->axon1 >= 0 && n->axon1 < cr->node_count ? cr->nodes[n->axon1]->value : 0.;
	float av2 = n->axon2 >= 0 && n->axon2 < cr->node_count ? cr->nodes[n->axon2]->value : 0.;
	if (n->op == 0) { // constant
	} else if (n->op == 1) { // time
		n->valueToBe = (float) sim->timer / 60.0;
	} else if (n->op == 2) { // x - coordinate
		n->valueToBe = n->x * 0.2;
	} else if (n->op == 3) { // y - coordinate
		n->valueToBe = -n->y * 0.2;
	} else if (n->op == 4) { // plus
		n->valueToBe = av1 + av2;
	} else if (n->op == 5) { // minus
		n->valueToBe = av1 - av2;
	} else if (n->op == 6) { // times
		n->valueToBe = av1 * av2;
	} else if (n->op == 7) { // divide
		n->valueToBe = av2 == 0. ? 0. : av1 / av2;
	} else if (n->op == 8) { // modulus
		n->valueToBe = av2 == 0. ? 0. : fmod(av1, av2);
	} else if (n->op == 9) { // sin
		n->valueToBe = sin(av1);
	} else if (n->op == 10) { // sig
		n->valueToBe = 1 / (1 + pow(2.71828182846, -av1));
	} else if (n->op == 11) { // pressure
		n->valueToBe = n->pressure;
	}
}

void node_realizeMathValues(struct node* n) {
	n->value = n->valueToBe;
}

struct node* node_copy(struct node* n) {
	return node_new(n->x, n->y, 0., 0., n->m, n->f, n->value, n->op, n->axon1, n->axon2);
}

void node_free(struct node* n) {
	free(n);
}

float evol_min(float f1, float f2) {
	return f1 > f2 ? f2 : f1;
}

float evol_max(float f1, float f2) {
	return f1 < f2 ? f2 : f1;
}

struct node* node_mutate(struct node* n, float mutability, struct creature* cr) {
	float newX = n->x + evol_r() * 0.5 * mutability;
	float newY = n->y + evol_r() * 0.5 * mutability;
	float newM = n->m + evol_r() * 0.1 * mutability;
	newM = evol_min(evol_max(newM, 0.3), 0.5);
	newM = 0.4;
	float newV = n->value * (1 + evol_r() * 0.2 * mutability);
	int newOperation = n->op;
	int newAxon1 = n->axon1;
	int newAxon2 = n->axon2;
	if (drand48() < bigMutationChance * mutability) {
		newOperation = (int) (drand48() * 12);
		if (newOperation == 12) newOperation = 11;
	}
	if (drand48() < bigMutationChance * mutability) {
		newAxon1 = (int) (drand48() * (float) cr->node_count);
		if (newAxon1 == cr->node_count) newAxon1--;
	}
	if (drand48() < bigMutationChance * mutability) {
		newAxon2 = (int) (drand48() * (float) cr->node_count);
		if (newAxon2 == cr->node_count) newAxon1--;
	}
	if (newOperation == 1) { // time
		newV = 0;
	} else if (newOperation == 2) { // x - coordinate
		newV = newX * 0.2;
	} else if (newOperation == 3) { // y - coordinate
		newV = -newY * 0.2;
	}
	return node_new(newX, newY, 0, 0, newM, evol_min(evol_max(n->f + evol_r() * 0.1 * mutability, 0), 1), newV, newOperation, newAxon1, newAxon2);
}

struct muscle* muscle_new(int axon, int c1, int c2, float len, float rigidity) {
	struct muscle* m = smalloc(sizeof(struct muscle));
	m->axon = axon;
	m->c1 = c1;
	m->c2 = c2;
	m->len = len;
	m->rigidity = rigidity;
	m->previousTarget = len;
	return m;
}

void muscle_applyForce(struct muscle* m, struct creature* cr, struct simulation* sim) {
	float target = m->previousTarget;
	if (energyDirection == 1 || (sim == NULL || sim->energy >= 0.0001)) {
		if (m->axon >= 0 && m->axon < cr->node_count) {
			target = cr->nodes[m->axon]->value;
			if (target < .5) target = .5;
			if (target > 1.5) target = 1.5;
			target *= m->len;
		} else {
			target = m->len;
		}
	}
	struct node* n1 = cr->nodes[m->c1];
	struct node* n2 = cr->nodes[m->c2];
	float dx = n2->x - n1->x;
	float dy = n2->y - n1->y;
	float distance = sqrt(dx * dx + dy * dy);
	float angle = atan2(n1->y - n2->y, n1->x - n2->x);
	float force = evol_min(evol_max(1. - (distance / target), -0.4), 0.4);
	n1->vx += cos(angle) * force * m->rigidity / n1->m;
	n1->vy += sin(angle) * force * m->rigidity / n1->m;
	n2->vx -= cos(angle) * force * m->rigidity / n2->m;
	n2->vy -= sin(angle) * force * m->rigidity / n2->m;
	if (sim != NULL) sim->energy = evol_max(sim->energy + energyDirection * abs(m->previousTarget - target) * m->rigidity * energyUnit, 0.);
	m->previousTarget = target;
}

struct muscle* muscle_copy(struct muscle* m) {
	return muscle_new(m->axon, m->c1, m->c2, m->len, m->rigidity);
}

struct muscle* muscle_mutate(struct muscle* m, float mutability, struct creature* cr) {
	int newc1 = m->c1;
	int newc2 = m->c2;
	int newAxon = m->axon;
	if (drand48() < bigMutationChance * mutability) {
		newc1 = (int) (drand48() * cr->node_count);
		if (newc1 == cr->node_count) newc1--;
	}
	if (drand48() < bigMutationChance * mutability) {
		newc2 = (int) (drand48() * cr->node_count);
		if (newc2 == cr->node_count) newc2--;
	}
	if (drand48() < bigMutationChance * mutability) {
		if (drand48() < .5) {
			newAxon = (int) (drand48() * cr->node_count);
			if (newAxon == cr->node_count) newAxon--;
		} else newAxon = -1;
	}
	float newR = evol_min(evol_max(m->rigidity * (1 + evol_r() * 0.9 * mutability), 0.01), 0.08);
	float newLen = evol_min(evol_max(m->len + evol_r() * mutability, 0.4), 1.25);
	return muscle_new(newAxon, newc1, newc2, newLen, newR);
}

struct creature* creature_new(int id, struct node** nodes, size_t node_count, struct muscle** muscles, size_t muscle_count, float d, int alive, float ct, float mut) {
	struct creature* c = smalloc(sizeof(struct creature));
	c->id = id;
	c->nodes = nodes;
	c->node_count = node_count;
	c->muscles = muscles;
	c->muscle_count = muscle_count;
	c->d = d;
	c->alive = alive;
	c->creatureTimer = ct;
	c->mutability = mut;
	return c;
}

void creature_addRandomMuscle(struct creature* cr, int c1, int c2) {
	int axon = -1;
	if (drand48() < .5) {
		axon = (int) (drand48() * (float) cr->node_count);
		if (axon == cr->node_count) axon--;
	}
	if (c1 == -1) {
		c1 = (int) (drand48() * (float) cr->node_count);
		if (c1 == cr->node_count) c1--;
		c2 = c1;
		while (c2 == c1 && cr->node_count >= 2) {
			c2 = (int) (drand48() * (float) cr->node_count);
			if (c2 == cr->node_count) c2--;
		}
	}
	float len = drand48() + .5;
	if (c1 != -1) {
		struct node* tc1 = cr->nodes[c1];
		struct node* tc2 = cr->nodes[c2];
		double dx = tc1->x - tc2->x;
		double dy = tc1->y - tc2->y;
		len = sqrt(dx * dx + dy * dy);
	}
	cr->muscles = srealloc(cr->muscles, (cr->muscle_count + 1) * sizeof(struct muscle*));
	cr->muscles[cr->muscle_count++] = muscle_new(axon, c1, c2, len, (drand48() * .06) + .02);
}

const int isim_operationAxons[] = { 0, 0, 0, 0, 2, 2, 2, 2, 2, 1, 1, 0 };

void creature_check(struct creature* nc) {
	//checkForOverlap();
	size_t nms = nc->muscle_count;
	for (int i = 0; i < nc->muscle_count; i++) {
		for (int j = i + 1; j < nc->muscle_count; j++) {
			struct muscle* mi = nc->muscles[i];
			struct muscle* mj = nc->muscles[j];
			if (mi->c1 == mj->c1 && mi->c2 == mj->c2) {
				free(mi);
				nc->muscles[i] = NULL;
				nms--;
				goto mcont;
			} else if (mi->c1 == mj->c2 && mi->c2 == mj->c1) {
				free(mi);
				nc->muscles[i] = NULL;
				nms--;
				goto mcont;
			} else if (mi->c1 == mi->c2 || mi->c1 == -1 || mi->c2 == -1) {
				free(mi);
				nc->muscles[i] = NULL;
				nms--;
				goto mcont;
			}
		}
		mcont: ;
	}
	if (nc->muscle_count > 0 && nc->muscles[nc->muscle_count - 1] != NULL && (nc->muscles[nc->muscle_count - 1]->c1 == -1 || nc->muscles[nc->muscle_count - 1]->c2 == -1)) {
		free(nc->muscles[nc->muscle_count - 1]);
		nc->muscles[nc->muscle_count - 1] = NULL;
		nms--;
	}
	if (nms != nc->muscle_count) {
		struct muscle** nm = smalloc(nms * sizeof(struct muscle*));
		size_t nmi = 0;
		for (size_t i = 0; i < nc->muscle_count; i++)
			if (nc->muscles[i] != NULL) nm[nmi++] = nc->muscles[i];
		free(nc->muscles);
		nc->muscles = nm;
		nc->muscle_count = nms;
	}
	//checkForLoneNodes();
	if (nc->node_count >= 3) {
		for (size_t i = 0; i < nc->node_count; i++) {
			int connections = 0;
			int connectedTo = -1;
			for (size_t j = 0; j < nc->muscle_count; j++) {
				if (nc->muscles[j]->c1 == i || nc->muscles[j]->c2 == i) {
					connections++;
					connectedTo = j;
				}
			}
			if (connections <= 1) {
				size_t newConnectionNode = (size_t)(drand48() * nc->node_count);
				if (newConnectionNode == nc->node_count) newConnectionNode--;
				while (newConnectionNode == i || newConnectionNode == connectedTo) {
					newConnectionNode = (size_t)(drand48() * nc->node_count);
					if (newConnectionNode == nc->node_count) newConnectionNode--;
				}
				creature_addRandomMuscle(nc, i, newConnectionNode);
			}
		}
	}
	//checkForBadAxons();
	for (size_t i = 0; i < nc->node_count; i++) { // this is probably unnecessary now
		struct node* ni = nc->nodes[i];
		if (ni->axon1 >= nc->node_count) {
			ni->axon1 = (size_t)(drand48() * nc->node_count);
			if (ni->axon1 == nc->node_count) ni->axon1--;
		}
		if (ni->axon2 >= nc->node_count) {
			ni->axon2 = (size_t)(drand48() * nc->node_count);
			if (ni->axon2 == nc->node_count) ni->axon2--;
		}
	}
	for (size_t i = 0; i < nc->muscle_count; i++) {
		struct muscle* mi = nc->muscles[i];
		if (mi->axon >= nc->node_count) {
			mi->axon = (size_t)(drand48() * nc->node_count);
			if (mi->axon == nc->node_count) mi->axon--;
		}
	}
	for (size_t i = 0; i < nc->node_count; i++) {
		struct node* ni = nc->nodes[i];
		ni->safeInput = (isim_operationAxons[ni->op] == 0);
	}
	int iterations = 0;
	int didSomething = 0;
	while (iterations < 1000) {
		didSomething = 0;
		for (size_t i = 0; i < nc->node_count; i++) {
			struct node* ni = nc->nodes[i];
			if (!ni->safeInput) {
				if ((isim_operationAxons[ni->op] == 1 && ni->axon1 >= 0 && ni->axon1 < nc->node_count && nc->nodes[ni->axon1]->safeInput) || (isim_operationAxons[ni->op] == 2 && ni->axon1 >= 0 && ni->axon1 < nc->node_count && nc->nodes[ni->axon1]->safeInput && ni->axon2 >= 0 && ni->axon2 < nc->node_count && nc->nodes[ni->axon2]->safeInput)) {
					ni->safeInput = 1;
					didSomething = 1;
				}
			}
		}
		if (!didSomething) iterations = 10000;
	}
	for (size_t i = 0; i < nc->node_count; i++) {
		struct node* ni = nc->nodes[i];
		if (!ni->safeInput) {
			ni->op = 0;
			ni->value = drand48();
		}
	}
}

struct creature* creature_copy(int id, struct creature* cr) {
	struct creature* cr2 = creature_new(id, smalloc(sizeof(struct node*) * cr->node_count), cr->node_count, smalloc(sizeof(struct muscle*) * cr->muscle_count), cr->muscle_count, cr->d, cr->alive, cr->creatureTimer, cr->mutability);
	for (size_t i = 0; i < cr->node_count; i++) {
		cr2->nodes[i] = node_copy(cr->nodes[i]);
	}
	for (size_t i = 0; i < cr->muscle_count; i++) {
		cr2->muscles[i] = muscle_copy(cr->muscles[i]);
	}
	return cr2;
}

struct creature* creature_modify(int id, struct creature* c) {
	struct creature* nc = creature_new(id, NULL, 0, NULL, 0, 0., 1, c->creatureTimer + evol_r() * 16. * c->mutability, evol_min(c->mutability * ((drand48() * .45) + .8), 2.));
	// i slightly changed the base game algorithm here, preselecting the node&muscle to be removed to optimize.
	ssize_t non = -1;
	ssize_t nom = -1;
	if (drand48() < bigMutationChance * c->mutability && c->node_count >= 4) {
		non = (ssize_t)(drand48() * c->node_count);
		if (non == c->node_count) non--;
	}
	if (drand48() < bigMutationChance * c->mutability && c->node_count >= 2) {
		nom = (ssize_t)(drand48() * c->muscle_count);
		if (nom == c->muscle_count) nom--;
	}
	nc->node_count = c->node_count - (non >= 0 ? 1 : 0);
	nc->nodes = smalloc(sizeof(struct node*) * c->node_count);
	nc->muscle_count = c->node_count - (nom >= 0 ? 1 : 0);
	nc->muscles = smalloc(sizeof(struct muscle*) * c->muscle_count);
	size_t ni = 0;
	for (size_t i = 0; i < c->node_count; i++) {
		if (i == non) continue;
		struct node* nn = node_mutate(c->nodes[i], c->mutability, c);
		if (non >= 0) {
			if (nn->axon1 >= non) nn->axon1--;
			if (nn->axon2 >= non) nn->axon2--;
		}
		nc->nodes[ni++] = nn;
	}
	ni = 0;
	for (size_t i = 0; i < c->muscle_count; i++) {
		if (i == nom) continue;
		struct muscle* nm = muscle_mutate(c->muscles[i], c->mutability, c);
		if (non >= 0) {
			if (nm->axon >= non) nm->axon--;
			if (nm->c1 >= non) nm->c1--;
			if (nm->c2 >= non) nm->c2--;
		}
		nc->muscles[ni++] = nm;
	}
	if (drand48() < bigMutationChance * c->mutability || c->node_count <= 2) {
		size_t pn = (size_t)(drand48() * nc->node_count);
		if (pn == nc->node_count) pn--;
		struct node* spn = nc->nodes[pn];
		float ang1 = drand48() * 2 * M_PI;
		float distance = sqrt(drand48());
		float x = spn->x + cos(ang1) * 0.5 * distance;
		float y = spn->y + sin(ang1) * 0.5 * distance;
		nc->nodes = srealloc(nc->nodes, (nc->node_count + 1) * sizeof(struct node*));
		int doc = (int) (drand48() * 12);
		if (doc == 12) doc--;
		int doc2 = (int) (drand48() * (nc->node_count + 1));
		if (doc2 == (nc->node_count + 1)) doc2--;
		int doc3 = (int) (drand48() * (nc->node_count + 1));
		if (doc3 == (nc->node_count + 1)) doc3--;
		nc->nodes[nc->node_count++] = node_new(x, y, 0, 0, 0.4, drand48(), drand48(), doc, doc2, doc3);
		int nextClosestNode = 0;
		float record = 100000;
		for (size_t i = 0; i < nc->node_count - 1; i++) {
			if (i != pn) {
				float dx = nc->nodes[i]->x - x;
				float dy = nc->nodes[i]->y - y;
				float d = sqrt(dx * dx + dy * dy);
				if (d < record) {
					record = d;
					nextClosestNode = i;
				}
			}
		}
		creature_addRandomMuscle(nc, pn, nc->node_count - 1);
		creature_addRandomMuscle(nc, nextClosestNode, nc->node_count - 1); // might be nice to make this a 90/10 or so
	}
	if (drand48() < bigMutationChance * c->mutability) {
		creature_addRandomMuscle(nc, -1, -1);
	}
	creature_check(nc);
	return nc;
}

void toStableConfiguration(struct creature* cr) {
	for (int j = 0; j < 200; j++) {
		for (size_t i = 0; i < cr->muscle_count; i++) {
			muscle_applyForce(cr->muscles[i], cr, NULL);
		}
		for (size_t i = 0; i < cr->node_count; i++) {
			node_applyForce(cr->nodes[i], NULL);
		}
	}
	for (size_t i = 0; i < cr->node_count; i++) {
		struct node* ni = cr->nodes[i];
		ni->vx = 0;
		ni->vy = 0;
	}
}

void adjustToCenter(struct creature* cr) {
	float avx = 0;
	float lowY = -1000;
	for (int i = 0; i < cr->node_count; i++) {
		struct node* ni = cr->nodes[i];
		avx += ni->x;
		if (ni->y + ni->m / 2. > lowY) {
			lowY = ni->y + ni->m / 2.;
		}
	}
	avx /= (float) cr->node_count;
	for (int i = 0; i < cr->node_count; i++) {
		struct node* ni = cr->nodes[i];
		ni->x -= avx;
		ni->y -= lowY;
	}
}

void mkcr_render(float partialTick) {
	glViewport(0, 0, width, height);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0., width, height, 0., 1000., 3000.);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0., 0., -2000.);
	glClearColor(1., 200. / 255., 130. / 255., 1.);
	char pbuf[256];
	snprintf(pbuf, 256, "Generation %i", (mkcr_gensel < 0 ? 0 : mkcr_gensel));
	text_left(pbuf, 20, 100, 6.);
	if (mkcr_gensel == -1) {
		glColor4f(100. / 255., 200. / 255., 100. / 255., 1.);
		rect(20, 250, 200, 100);
		glColor4f(1., 1., 1., 1.);
		text_left("Since there are no creatures yet, create 1000 creatures!", 20, 160, 4.);
		text_left("They will be randomly created, and also very simple.", 20, 200, 4.);
		text_center("CREATE", 120, 325, 6.);
	} else {
		glColor4f(100. / 255., 200. / 255., 100. / 255., 1.);
		rect(760, 20, 460, 40);
		rect(760, 70, 460, 40);
		rect(760, 120, 230, 40);
		//if(gensToDo >= 2) {
		// 		glColor4f(128. / 255., 200. / 255., 100 / 128., 1.);
		//}else{
		glColor4f(70. / 255., 140. / 255., 70. / 255., 1.);
		//}
		rect(990, 120, 230, 40);
		text_left("Do 1 step-by-step generation.", 770, 50, 4.);
		text_left("Do 1 quick generation.", 770, 100, 4.);
		text_left("Do 1 gen ASAP.", 770, 150, 4.);
		text_left("Do gens ALAP.", 1000, 150, 4.);
		snprintf(pbuf, 256, "Median %s", fitnessName);
		text_left(pbuf, 50, 160, 4.);
		//text_right(float(round(percentile.get(min(genSelected, percentile.size()-1))[14]*1000))/1000+" "+fitnessUnit, 700, 160);
		drawHistogram(&mkcr_mhg, 760, 410, 460, 280);
		drawGraphImage(&mkcr_mg, 70, 180, 650, 380, 100., -10.);
		drawSegBars(70, 580, 650, 100);
		if (mkcr_gensel >= 1) {
			//float genWidth = 590.0 / (float) (mkcr_gen == 0 ? 1 : mkcr_gen);
			//float lineX = 110. + (float) mkcr_gensel * genWidth;
			//line(lineX, 180, lineX, 500 + 180);
			/*Integer[] s = speciesCounts.get(genSelected);
			 //textAlign (LEFT);
			 //textFont(font, 12);
			 noStroke();
			 for (int i = 1; i < 101; i++) {
			 int c = s[i] - s[i - 1];
			 if (c >= 25) {
			 float y = ((s[i] + s[i - 1]) / 2) / 1000.0 * 100 + 573;
			 if (i - 1 == topSpeciesCounts.get(genSelected)) {
			 stroke(0);
			 strokeWeight(2);
			 } else {
			 noStroke();
			 }
			 fill(255, 255, 255);
			 rect(lineX + 3, y, 56, 14);
			 colorMode(HSB, 1.0);
			 fill(getColor(i - 1, true));
			 text("S" + floor((i - 1) / 10) + "" + ((i - 1) % 10) + ": " + c, lineX + 5, y + 11);
			 colorMode(RGB, 255);
			 }
			 }
			 noStroke();*/
		}
	}
}

struct simulation isim;

float isim_averageX;
float isim_averageY;

void isim_setAverages(struct creature* cr) {
	isim_averageX = 0;
	isim_averageY = 0;
	for (size_t i = 0; i < cr->node_count; i++) {
		struct node* ni = cr->nodes[i];
		isim_averageX += ni->x;
		isim_averageY += ni->y;
	}
	isim_averageX = isim_averageX / (float) cr->node_count;
	isim_averageY = isim_averageY / (float) cr->node_count;
}

void creature_simulate(struct creature* cr, struct simulation* sim) {
	for (size_t i = 0; i < cr->muscle_count; i++) {
		muscle_applyForce(cr->muscles[i], cr, sim);
	}
	for (size_t i = 0; i < cr->node_count; i++) {
		struct node* ni = cr->nodes[i];
		node_applyGravity(ni);
		node_applyForce(ni, sim);
		node_hitWalls(ni, sim);
		node_doMath(ni, sim, cr);
	}
	for (size_t i = 0; i < cr->node_count; i++) {
		node_realizeMathValues(cr->nodes[i]);
	}
	sim->averageNodeNausea = sim->totalNodeNausea / (float) cr->node_count;
	sim->timer++;
}

float isim_lf[1000];

void mkcr_button(int action, int mods, double x, double y) {
	if (action == 0) {
		for (int y = 0; y < 25; y++) {
			for (int x = 0; x < 40; x++) {
				int nodeNum = (int) (drand48() * 4.) + 3;
				int muscleNum = (int) (drand48() * (float) (nodeNum * 2 - 5)) + nodeNum + 1;
				struct node** nodes = smalloc(nodeNum * sizeof(struct node*));
				struct muscle** muscles = smalloc(muscleNum * sizeof(struct muscle*));
				for (int i = 0; i < nodeNum; i++) {
					nodes[i] = node_new(drand48() * 2. - 1., drand48() * 2. - 1., 0, 0, 0.4, drand48(), drand48(), (int) (drand48() * 12), (int) (drand48() * nodeNum), (int) (drand48() * nodeNum));
				}
				for (int i = 0; i < muscleNum; i++) {
					int tc1 = 0;
					int tc2 = 0;
					int taxon = -1;
					if (drand48() < .5) taxon = (int) (drand48() * nodeNum);
					if (i < nodeNum - 1) {
						tc1 = i;
						tc2 = i + 1;
					} else {
						tc1 = (int) (drand48() * nodeNum);
						tc2 = tc1;
						while (tc2 == tc1) {
							tc2 = (int) (drand48() * nodeNum);
						}
					}
					float s = 0.8;
					if (i >= 10) {
						s *= 1.414;
					}
					float len = drand48() + .5;
					muscles[i] = muscle_new(taxon, tc1, tc2, len, drand48() * .06 + .02);
				}
				float heartbeat = drand48() * 40. + 40.;
				cpop[y * 40 + x] = creature_new(y * 40 + x + 1, nodes, nodeNum, muscles, muscleNum, 0, 1, heartbeat, 1.0);
				toStableConfiguration(cpop[y * 40 + x]);
				adjustToCenter(cpop[y * 40 + x]);
//drawCreature(cpop[y * 40 + x], x * 3 + 5.5, y * 2.5 + 3, 0);
				creature_check(cpop[y * 40 + x]);
			}
		}
		guistate_set(2);
	} else if (action == 1 || action == 2 || action == 3 || action == 4 || action == 5) {
		for (size_t i = 0; i < 1000; i++) {
			ccsimc[i] = creature_copy(cpop[i]->id, cpop[i]);
		}
		isim.averageNodeNausea = 0.;
		isim.energy = 0.;
		isim.timer = 0;
		isim.totalNodeNausea = 0.;
		if (action == 1) {
			guistate_set(3);
			return;
		} else if (action == 2 || action == 3) {
			for (int i = 0; i < 1000; i++) {
				for (int j = 0; j <= 900; j++) {
					creature_simulate(ccsimc[i], &isim);
					isim_setAverages(ccsimc[i]);
				}
				isim_lf[i] = -isim_averageX;
				isim.timer = 0;
				isim.energy = 0.;
				isim.averageNodeNausea = 0.;
				isim.totalNodeNausea = 0.;
			}
			if (action == 2) guistate_set(2);
			else {
				struct pqueue* sq = new_pqueue(1000, 0);
				for (int i = 0; i < 1000; i++) {
					add_pqueue(sq, cpop[i], isim_lf[i]);
				}
				for (int i = 0; i < 1000; i++) {
					cpop[i] = pop_pqueue(sq);
				}
				del_pqueue(sq);
				//
				for (int j = 0; j < 500; j++) {
					float f = (float) j / 1000.;
					float rand = (pow((drand48() * 2. - 1.), 3.) + 1.) / 2.; //cube function
					int slowDies = (f <= rand);
					int j2;
					int j3;
					if (slowDies) {
						j2 = j;
						j3 = 999 - j;
					} else {
						j2 = 999 - j;
						j3 = j;
					}
					struct creature* cj = cpop[j2];
					cj->alive = 1;
					struct creature* ck = cpop[j3];
					ck->alive = 0;
				}
				//
				for (int j = 0; j < 500; j++) {
					int j2 = j;
					if (!cpop[j]->alive) j2 = 999 - j;
					struct creature* cj = cpop[j2];
					struct creature* cj2 = cpop[999 - j2];
					cpop[j2] = creature_copy(cj->id + 1000, cj);
					cpop[999 - j2] = creature_modify(cj2->id + 1000, cj);
					toStableConfiguration(cpop[999 - j2]);
					adjustToCenter(cpop[999 - j2]);
				}
				mkcr_gen++;
				mkcr_gensel = mkcr_gen;
			}
		}
	}
}

void mkcr_init() {
	mkcr_mhg.data = NULL;
	mkcr_mhg.data_count = 0;
	mkcr_mhg.median = 0.;
	mkcr_mg.data = NULL;
	mkcr_mg.data_count = 0;
	mkcr_mg.data_subcount = 0;
	setFont (TX_FONT);
	struct gui_button btn;
	btn.action = 0;
	btn.width = 200;
	btn.height = 100;
	btn.x = 20;
	btn.y = 250;
	gui_addbutton(1, btn);
}

float isim_camX;
float isim_camY;
float isim_camZoom;

void isim_drawGround(struct simulation* sim) {
	int stairDrawStart = evol_max(1, (int) (-isim_averageY / hazelStairs) - 10);
	glColor4f(0., 130. / 255., 0., 1.);
	if (haveGround) {
		//rect(isim_camX - width / 2., isim_camY, width, height / 2.);
		rect(isim_camX - width / 2., 0, width, height);
	}
	for (size_t i = 0; i < sim->rect_count; i++) {
		struct sim_rect sr = sim->rects[i];
		rect(sr.x1, sr.y1, sr.x2 - sr.x1, sr.y2 - sr.y1);
	}
	if (hazelStairs > 0) {
		for (int i = stairDrawStart; i < stairDrawStart + 20; i++) {
			glColor4f(1., 1., 1., .5);
			rect((isim_averageX - 20), -hazelStairs * i, 40, hazelStairs * 0.3);
			glColor4f(1., 1., 1., 1.);
			rect((isim_averageX - 20), -hazelStairs * i, 40, hazelStairs * 0.15);
		}
	}
}

void drawNode(struct node* ni, float x, float y) {
	if (ni->f <= 0.5) glColor4f(1., (255. - ni->f * 512.) / 255., (255. - ni->f * 512.) / 255., 1.);
	else glColor4f((512. - ni->f * 512.) / 255., 0., 0., 1.);
	//ellipse((ni.x + x), (ni.y + y), ni.m, ni.m);
	circle(ni->x + x, ni->y + y, ni->m);
	if (ni->f >= 0.5) {
		glColor4f(1., 1., 1., 1.);
	} else {
		glColor4f(0., 0., 0., 1.);
	}
	//rect(ni->x + x - ni->m / 2., ni->y + y - ni->m / 2., ni->m / 2., ni->m / 2.);
	static char* operationNames[] = { "#", "time", "px", "py", "+", "-", "*", "÷", "%", "sin", "sig", "pres" };
	//char pbuf[256];
	//snprintf(pbuf, 256, "%.2f", ni->value);
	//text_center(pbuf, (ni->x + x + ni->m / 2.), (ni->y + ni->m / 2. * lineY2 + y), .03);
	//text_center(operationNames[ni->op], (ni->x + x + ni->m / 2.), (ni->y + ni->m / 2. * lineY1 + y), .03);
}

void drawSingleAxon(float x1, float y1, float x2, float y2) {
	float arrowHeadSize = 0.1;
	float angle = atan2(y2 - y1, x2 - x1);
	glColor4f(1., 1., 0., 1.);
	line(x1, y1, x2, y2);
	line(x1, y1, (x1 + cos(angle + M_PI * 0.25) * arrowHeadSize), (y1 + sin(angle + M_PI * 0.25) * arrowHeadSize));
	line(x1, y1, (x1 + cos(angle + M_PI * 1.75) * arrowHeadSize), (y1 + sin(angle + M_PI * 1.75) * arrowHeadSize));
}

void drawNodeAxons(struct node** nodes, size_t node_count, int i, float x, float y) {
	struct node* ni = nodes[i];
	if (isim_operationAxons[ni->op] >= 1) {
		struct node* axonSource = nodes[ni->axon1];
		float point1x = ni->x - ni->m * 0.3 + x;
		float point1y = ni->y - ni->m * 0.3 + y;
		float point2x = axonSource->x + x;
		float point2y = axonSource->y + axonSource->m * 0.5 + y;
		drawSingleAxon(point1x, point1y, point2x, point2y);
	}
	if (isim_operationAxons[ni->op] == 2) {
		struct node* axonSource = nodes[ni->axon1];
		float point1x = ni->x + ni->m * 0.3 + x;
		float point1y = ni->y - ni->m * 0.3 + y;
		float point2x = axonSource->x + x;
		float point2y = axonSource->y + axonSource->m * 0.5 + y;
		drawSingleAxon(point1x, point1y, point2x, point2y);
	}
}

void drawMuscle(struct muscle* mi, struct node** nodes, size_t node_count, float x, float y) {
	struct node* ni1 = nodes[mi->c1];
	struct node* ni2 = nodes[mi->c2];
	float w = 0.15;
	if (mi->axon >= 0 && mi->axon < node_count) {
		w = nodes[mi->axon]->value;
		if (w > 1.5) w = 1.5;
		if (w < .5) w = .5;
		w *= 0.15;
	}
	glColor4f(70. / 255., 35. / 255., 0., (mi->rigidity / 3000.) / 255.);
	line((ni1->x + x), (ni1->y + y), (ni2->x + x), (ni2->y + y));
}

void drawMuscleAxons(struct muscle* mi, struct node** nodes, size_t node_count, float x, float y) {
	struct node* ni1 = nodes[mi->c1];
	struct node* ni2 = nodes[mi->c2];
	if (mi->axon >= 0 && mi->axon < node_count) {
		struct node* axonSource = nodes[mi->axon];
		float muscleMidX = (ni1->x + ni2->x) * 0.5 + x;
		float muscleMidY = (ni1->y + ni2->y) * 0.5 + y;
		drawSingleAxon(muscleMidX, muscleMidY, axonSource->x + x, axonSource->y + axonSource->m * 0.5 + y);
		float averageMass = (ni1->m + ni2->m) * 0.5;
		glColor4f(1., 1., 0., 1.);
		char pbuf[256];
		float w = nodes[mi->axon]->value;
		if (w > 1.5) w = 1.5;
		if (w < .5) w = .5;
		snprintf(pbuf, 256, "%.2f", w);
		//text_center(pbuf, muscleMidX, muscleMidY, 2.); // todo incorporate averageMass into scale
	}
}

void drawPosts() {
	int startPostY = evol_min(-8, (int) (isim_averageY / 4) * 4 - 4);
	for (int postY = startPostY; postY <= startPostY + 8; postY += 4) {
		for (int i = (int) (isim_averageX / 5 - 5); i <= (int) (isim_averageX / 5 + 5); i++) {
			glColor4f(1., 1., 1., 1.);
			rect((i * 5.0 - 0.1), (-3.0 + postY), 0.2, 3.0);
			rect((i * 5.0 - 1), (-3.0 + postY), 2.0, 1.0);
			glColor4f(120. / 255., 120. / 255., 120. / 255., 1.);
			char pbuf[256];
			snprintf(pbuf, 256, "%i m", i);
			//text_center(pbuf, i * 5.0, (-2.17 + postY), 4.);
		}
	}
}

void drawCreature(struct creature* cj, float x, float y) {
	for (size_t i = 0; i < cj->muscle_count; i++) {
		drawMuscle(cj->muscles[i], cj->nodes, cj->node_count, x, y);
	}
	for (size_t i = 0; i < cj->node_count; i++) {
		drawNode(cj->nodes[i], x, y);
	}
	for (size_t i = 0; i < cj->muscle_count; i++) {
		drawMuscleAxons(cj->muscles[i], cj->nodes, cj->node_count, x, y);
	}
	for (size_t i = 0; i < cj->node_count; i++) {
		drawNodeAxons(cj->nodes, cj->node_count, i, x, y);
	}
}

int vcrs_state = 0;

struct vcrs_tcreat {
		struct creature* cr;
		uint16_t i;
};

struct vcrs_tcreat vcrs_npi[1000];
float vcrs_interpcr_x[1000];
float vcrs_interpcr_y[1000];
int vcrs_tick = 60;

void vcrs_render(float partialTick) {
	glViewport(0, 0, width, height);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0., width, height, 0., 1000., 3000.);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0., 0., -2000.);
	glClearColor(220. / 255., 253. / 255., 102. / 255., 1.);
	glPushMatrix();
	glScalef(10., 10., 10.);
	float pa = (float) vcrs_tick / 60.;
	for (int y = 0; y < 25; y++) {
		for (int x = 0; x < 40; x++) {
			if (cpop[y * 40 + x]->alive) {
				if (vcrs_tick < 60) {
					glPushMatrix();
					glTranslatef(vcrs_interpcr_x[y * 40 + x], vcrs_interpcr_y[y * 40 + x], 0.);
					glRotatef(pa * 360., 0., 0., 1.);
					drawCreature(cpop[y * 40 + x], 0., 0.);
					glPopMatrix();
				} else drawCreature(cpop[y * 40 + x], (float) x * 3. + 5.5, (float) y * 2.5 + 3.);
			} else {
				glColor4f(0., 0., 0., 1.);
				rect(x * 3. + 4.5, y * 2.5 + 2., 3., 2.5);
			}
		}
	}
	glPopMatrix();
	glColor4f(100. / 255., 100. / 255., 200. / 255., 1.);
	rect(900, 664, 260, 40);
	text_center("Here are your 1000 randomly generated creatures!", 1280 / 2 - 200, 700, 4.);
	if (mkcr_gensel == -1) {
		text_center("Back", 1280 - 250, 700, 4.);
	} else if (vcrs_state == 0) {
		text_center("Sort", 1280 - 250, 700, 4.);
	} else if (vcrs_state == 1) {
		text_center("Kill 500", 1280 - 250, 700, 4.);
	} else if (vcrs_state == 2) {
		text_center("Reproduce", 1280 - 250, 700, 4.);
	} else if (vcrs_state == 3) {
		text_center("Back", 1280 - 250, 700, 4.);
	}
}

void vcrs_button(int action, int mods, double x, double y) {
	if (action == 0) {
		if (mkcr_gensel == -1) {
			gui_rembutton(1, 0);
			struct gui_button bt;
			bt.action = 1;
			bt.x = 760;
			bt.y = 20;
			bt.width = 460;
			bt.height = 40;
			gui_addbutton(1, bt);
			bt.action = 2;
			bt.x = 760;
			bt.y = 70;
			bt.width = 460;
			bt.height = 40;
			gui_addbutton(1, bt);
			bt.action = 3;
			bt.x = 760;
			bt.y = 120;
			bt.width = 230;
			bt.height = 40;
			gui_addbutton(1, bt);
		} else if (vcrs_state == 0) {
			struct pqueue* sq = new_pqueue(1000, 0);
			for (int i = 0; i < 1000; i++) {
				struct vcrs_tcreat* vc = smalloc(sizeof(struct vcrs_tcreat));
				vc->cr = cpop[i];
				vc->i = i;
				add_pqueue(sq, vc, isim_lf[i]);
			}
			for (int i = 0; i < 1000; i++) {
				struct vcrs_tcreat* vt = pop_pqueue(sq);
				cpop[i] = vt->cr;
				vcrs_npi[i] = *vt;
				free(vt);
			}
			del_pqueue(sq);
			vcrs_state = 1;
			vcrs_tick = 0;
			return;
		} else if (vcrs_state == 1) {
			for (int j = 0; j < 500; j++) {
				float f = (float) j / 1000.;
				float rand = (pow((drand48() * 2. - 1.), 3.) + 1.) / 2.; //cube function
				int slowDies = (f <= rand);
				int j2;
				int j3;
				if (slowDies) {
					j2 = j;
					j3 = 999 - j;
				} else {
					j2 = 999 - j;
					j3 = j;
				}
				struct creature* cj = cpop[j2];
				cj->alive = 1;
				struct creature* ck = cpop[j3];
				ck->alive = 0;
			}
			vcrs_state = 2;
			return;
		} else if (vcrs_state == 2) {
			for (int j = 0; j < 500; j++) {
				int j2 = j;
				if (!cpop[j]->alive) j2 = 999 - j;
				struct creature* cj = cpop[j2];
				struct creature* cj2 = cpop[999 - j2];
				cpop[j2] = creature_copy(cj->id + 1000, cj);
				cpop[999 - j2] = creature_modify(cj2->id + 1000, cj);
				toStableConfiguration(cpop[999 - j2]);
				adjustToCenter(cpop[999 - j2]);
			}
			//for (int j = 0; j < 1000; j++) {
			//	Creature cj = c2.get(j);
			//	c[cj.id - (gen * 1000) - 1001] = cj.copyCreature(-1);
			//}
			mkcr_gen++;
			vcrs_state = 3;
			return;
		}
		vcrs_state = 0;
		mkcr_gensel++;
		mkcr_gensel = mkcr_gen;
		guistate_set(1);
	}
}

void vcrs_stick() {
	if (vcrs_tick < 60) {
		float pa = (float) vcrs_tick / 60.;
		for (int y = 0; y < 25; y++) {
			for (int x = 0; x < 40; x++) {
				int ni = y * 40 + x;
				int oi = vcrs_npi[ni].i;
				float npx = (float) x * 3. + 5.5;
				float npy = (float) y * 2.5 + 3.;
				float opx = (float) (oi % 40) * 3. + 5.5;
				float opy = (float) (oi / 40) * 2.5 + 3.;
				vcrs_interpcr_x[ni] = opx * (1. - pa) + npx * pa;
				vcrs_interpcr_y[ni] = opy * (1. - pa) + npy * pa;
			}
		}
		vcrs_tick++;
	}
}

void vcrs_init() {
	mkcr_mhg.data = NULL;
	mkcr_mhg.data_count = 0;
	mkcr_mhg.median = 0.;
	mkcr_mg.data = NULL;
	mkcr_mg.data_count = 0;
	mkcr_mg.data_subcount = 0;
	setFont (TX_FONT);
	struct gui_button btn;
	btn.action = 0;
	btn.width = 260;
	btn.height = 40;
	btn.x = 900;
	btn.y = 664;
	gui_addbutton(2, btn);
}

size_t ccsmc = 0;
int ccspd = 1;

void drawStats(float x, float y, float size) {
	glPushMatrix();
	glTranslatef(x, y, 0.);
	glScalef(size, size, size);
	char pbuf[256];
	snprintf(pbuf, 256, "Creature ID: %i", ccsimc[ccsmc]->id);
	text_right(pbuf, 0, 32, 4.);
	float time = 0.;
	if (ccspd > 60) {
		time = fmod((float) ((isim.timer + ccsmc * 37) / 60.), 15);
	} else {
		time = ((float) isim.timer / 60.);
	}
	snprintf(pbuf, 256, "Time: %.2f", time);
	text_right(pbuf, 0, 64, 4.);
	snprintf(pbuf, 256, "Playback Speed: %i", (int) evol_max(1., ccspd));
	text_right(pbuf, 0, 96, 4.);
	/*String extraWord = "used";
	 if (energyDirection == -1) {
	 extraWord = "left";
	 }*/
	snprintf(pbuf, 256, "X: %.2f", isim_averageX / 5.);
	text_right(pbuf, 0, 128, 4.);
	snprintf(pbuf, 256, "Y: %.2f", isim_averageY / 5.);
	text_right(pbuf, 0, 160, 4.);
	snprintf(pbuf, 256, "Energy %s: %.2f yums", (energyDirection == -1 ? "left" : "used"), isim.energy);
	text_right(pbuf, 0, 192, 4.);
	snprintf(pbuf, 256, "A.N.Nausea %.2f blehs", isim.averageNodeNausea);
	text_right(pbuf, 0, 224, 4.);
	glPopMatrix();
}

void isim_render(float partialTick) {
	glViewport(0, 0, width, height);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0., width, height, 0., 1000., 3000.);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0., 0., -2000.);
	glClearColor(120. / 255., 200. / 255., 255. / 255., 1.);
	glPushMatrix();
	glTranslatef((float) width / 2.0, (float) height / 2.0, 0.);
//glScalef(1.0 / isim_camZoom, 1.0 / isim_camZoom, 1.0 / isim_camZoom);
	glScalef(25., 25., 25.);
	glTranslatef(-isim_camX, -isim_camY, 0.);
	drawPosts();
	isim_drawGround(&isim);
	drawCreature(ccsimc[ccsmc], 0, 0);
//drawArrow (averageX);
	glPopMatrix();
	drawStats(width - 10, 0, 0.7);
	glColor4f(0., 0., 0., 1.);
	rect(0, 0, 90, 40);
	rect(120, 0, 120, 40);
//rect(120, height - 40, 240, 40);
//text_center_white
	glEnable (GL_TEXTURE_2D);
	glEnable (GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//
	glPushMatrix();
	glTranslatef(45., 35., 0.);
	glScalef(4., 4., 4.);
	drawString("SKIP", -((float) stringWidth("SKIP") / 2.), -10., 0xFFFFFFFF);
	glPopMatrix();
//
	glPushMatrix();
	glTranslatef(175., 35., 0.);
	glScalef(4., 4., 4.);
	drawString("FINISH", -((float) stringWidth("FINISH") / 2.), -10., 0xFFFFFFFF);
	glPopMatrix();
//
	glDisable(GL_BLEND);
	glDisable(GL_TEXTURE_2D);

}

void isim_button(int action, int mods, double x, double y) {
	if (action == 0) {
		for (int i = isim.timer; i <= 900; i++) {
			creature_simulate(ccsimc[ccsmc], &isim);
		}
		isim_setAverages(ccsimc[ccsmc]);
		isim_lf[ccsmc] = -isim_averageX;
		ccsmc++;
		if (ccsmc >= 1000) {
			ccsmc = 0;
			guistate_set(2);
		}
		isim_camX = 0.;
		if (ccspd < 1000) ccspd++;
		isim.timer = 0;
		isim.energy = 0.;
		isim.averageNodeNausea = 0.;
		isim.totalNodeNausea = 0.;
	} else if (action == 1) {
		while (ccsmc < 1000) {
			for (int i = isim.timer; i <= 900; i++) {
				creature_simulate(ccsimc[ccsmc], &isim);
			}
			isim_setAverages(ccsimc[ccsmc]);
			isim_lf[ccsmc] = -isim_averageX;
			ccsmc++;
			if (ccsmc >= 1000) break;
			isim_camX = 0.;
			if (ccspd < 1000) ccspd++;
			isim.timer = 0;
			isim.energy = 0.;
			isim.averageNodeNausea = 0.;
			isim.totalNodeNausea = 0.;
		}
		ccspd = 1;
		ccsmc = 0;
		guistate_set(2);
	}
}

void isim_init() {
//set rects in isim here
	setFont (TX_FONT);
	struct gui_button btn;
	btn.action = 0;
	btn.x = 0;
	btn.y = 0;
	btn.width = 90;
	btn.height = 40;
	gui_addbutton(3, btn);
	btn.action = 1;
	btn.x = 120;
	btn.y = 0;
	btn.width = 120;
	btn.height = 40;
	gui_addbutton(3, btn);
}

void isim_tick() {
	for (int i = 0; i < ccspd; i++) {
		if (isim.timer < 900) {
			creature_simulate(ccsimc[ccsmc], &isim);
		}
	}
	isim_setAverages(ccsimc[ccsmc]);
	if (ccspd < 30) {
		for (int s = 0; s < ccspd; s++) {
			isim_camX += (isim_averageX - isim_camX) * 0.06;
			isim_camY += (isim_averageY - isim_camY) * 0.06;
		}
	} else {
		isim_camX = isim_averageX;
		isim_camY = isim_averageY;
	}
	if (isim.timer == 900) {
		isim_lf[ccsmc] = -isim_averageX;
		ccsmc++;
		if (ccsmc >= 1000) {
			ccspd = 1;
			ccsmc = 0;
			guistate_set(2);
		}
		isim_camX = 0.;
		if (ccspd < 1000) ccspd++;
		isim.timer = 0;
		isim.energy = 0.;
		isim.averageNodeNausea = 0.;
		isim.totalNodeNausea = 0.;
	}
}
