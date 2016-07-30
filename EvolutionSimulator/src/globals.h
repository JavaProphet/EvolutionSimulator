/*
 * config.h
 *
 *  Created on: Jul 25, 2016
 *      Author: root
 */

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <GLFW/glfw3.h>
#include <stdlib.h>

int tps;
int width;
int height;
GLFWwindow* window;
int frames;
char* windowTitle;
uint32_t frameLimit;
int swidth;
int sheight;
int csf; // current scale factor
char* fitnessName;
char* fitnessUnit;
int minBar;
int maxBar;
int barLen;
int histBarsPerMeter;
float airFriction;
float gravity;
float nauseaUnit;
float pressureUnit;
float friction;
float energyUnit;
int haveGround;
float hazelStairs;
float bigMutationChance;
int energyDirection;
float lineY1;
float lineY2;

#endif /* GLOBALS_H_ */
