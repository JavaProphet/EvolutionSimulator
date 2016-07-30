/*
 * evolsim.h
 *
 *  Created on: Jul 26, 2016
 *      Author: root
 */

#ifndef EVOLSIM_H_
#define EVOLSIM_H_

void evol_render(float partialTick);

void evol_init();

void evol_button(int action, int mods, double x, double y);

void mkcr_render(float partialTick);

void mkcr_button(int action, int mods, double x, double y);

void mkcr_init();

void vcrs_render(float partialTick);

void vcrs_button(int action, int mods, double x, double y);

void vcrs_init();

void vcrs_stick();

void isim_render(float partialTick);

void isim_button(int action, int mods, double x, double y);

void isim_init();

void isim_tick();

#define TX_FONT 1

#endif /* EVOLSIM_H_ */
