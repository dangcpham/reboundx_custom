/**
 * @file    core.h
 * @brief   Central internal functions for REBOUNDx (not called by user)
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef REBX_CORE_H
#define REBX_CORE_H

#include <stdint.h>
#include "rebound.h"

/****************************************
Basic types in REBOUNDx
*****************************************/

/* 	Main structure used for all parameters added to particles.
 	These get added as nodes to a linked list for each particle, stored at particles[i].ap.*/
struct rebx_param{
    void* paramPtr;                     // Pointer to the parameter (void* so it can point to different types).
    uint32_t hash;                      // Hash for the parameter name.
    struct rebx_param* next;            // Pointer to the next parameter in the linked list.
};

/*  Structure for all REBOUNDx effects.
 *  These get added as nodes to the effects linked list in the rebx_extras structure.*/
struct rebx_effect{
	rebx_param* ap;					    // Linked list of parameters for the effect.
    void (*force) (struct reb_simulation* sim, struct rebx_effect* effect); // Pointer to function to call during forces evaluation.
    void (*ptm) (struct reb_simulation* sim, struct rebx_effect* effect);   // Pointer to function to call after each timestep.
	struct rebx_effect* next;			// Pointer to the next effect in the linked list.
};

/*	Nodes for a linked list to all the parameters that have been allocated by REBOUNDx (so it can later free them).*/
struct rebx_param_to_be_freed{
    struct rebx_param* param;           // Pointer to a parameter node allocated by REBOUNDx.
    struct rebx_param_to_be_freed* next;// Pointer to the next node in the linked list rebx_extras.params_to_be_freed.
};

/****************************************
Main REBOUNDx structure
*****************************************/
struct rebx_extras {	
	struct reb_simulation* sim;								// Pointer to the simulation REBOUNDx is linked to.
	struct rebx_effect* effects;		                    // Linked list with pointers to all the effects added to the simulation.
	struct rebx_param_to_be_freed* params_to_be_freed; 		// Linked list with pointers to all parameters allocated by REBOUNDx (for later freeing).

};

/*****************************
 Internal initialization routine.
 ****************************/

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx); // Initializes all pointers and values.

/*****************************
 Garbage Collection Routines
 ****************************/

void rebx_free_params(struct rebx_extras* rebx);            // Steps through linked list to free all allocated particle parameters.
void rebx_free_effects(struct rebx_extras* rebx);           // Frees all effects in effects linked list 

/**********************************************
 Functions executing forces & ptm each timestep
 *********************************************/

void rebx_forces(struct reb_simulation* sim);                       // Calls all the forces that have been added to the simulation.
void rebx_post_timestep_modifications(struct reb_simulation* sim);  // Calls all the post-timestep modifications that have been added to the simulation.

/**********************************************
 Adders for linked lists in extras structure
 *********************************************/

void rebx_add_effect(struct rebx_extras* rebx, const char* name);

// Add a parameter to the params_to_be_freed linked list for later freeing.
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param); // add a node for param in the rebx_params_to_be_freed linked list.

/*********************************************************************************
 General particle parameter getter
 ********************************************************************************/

void* rebx_get_ap_hash(struct rebx_param current, uint32_t hash);   // Returns rebx_param corresponding to hash.  If it doesn't exist, returns NULL.

/*********************************************************************************
 Getters and Setters for particle parameters (need new set for each variable type)
 ********************************************************************************/
void rebx_set_ap_double_hash(struct rebx_extras* rebx, struct rebx_param* ap, uint32_t hash, double value);
double* rebx_add_ap_double(struct rebx_extras* rebx, struct rebx_param* ap, uint32_t hash);               
double rebx_get_ap_double_hash(struct rebx_param* ap, uint32_t hash);                                    

/***********************************************************************************
 * Miscellaneous Functions
***********************************************************************************/

double install_test(void);  // Function for testing whether REBOUNDx can load librebound.so and call REBOUND functions. */

#endif
