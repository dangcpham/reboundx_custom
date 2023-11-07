/** * @file planet_force.c
 * @brief   Force due to planets for massless particles
 * @author  Dang Pham
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
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Central Force$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 Dang Pham
 * Implementation Paper    
 * Based on                None
 * C Example               None
 * Python Example          None.
 * ======================= ===============================================
 * 
 * Adds the effect of a planet on an inclined circular orbit on massless particles.
 * 
 *
 *
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

void rebx_planets_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    struct rebx_planets* const planets = rebx_get_param(rebx, force->ap, "planets_data");
    const double Nplanets = planets->Nplanets;
    const double Nparams = planets->Nparams;
    const double t = sim->t;
    const double G = sim->G;

    int idx = rc_idx(1, 2, Nparams);
    printf("%e\n", planets->planets_data[idx]);
}

int rc_idx(const int r, const int c, const int N){
    return r*N + c;
}

struct rebx_planets* rebx_create_planets(struct rebx_extras* const rebx, const int Nplanets, const int Nparams, const double* planets_data, const double* masses ){
    struct rebx_planets* planets = rebx_malloc(rebx, sizeof(*planets));
    rebx_init_planets(rebx, planets, Nplanets, Nparams, planets_data, masses);
    return planets;
}

void rebx_init_planets(struct rebx_extras* const rebx, struct rebx_planets* const planets, const int Nplanets, const int Nparams, const double* planets_data, const double* masses){
    planets->Nplanets = Nplanets;
    planets->Nparams = Nparams;
    planets->planets_data = calloc(Nplanets * Nparams, sizeof(*planets->planets_data));
    planets->masses = calloc(Nplanets + 1, sizeof(*planets->masses));

    memcpy(planets->planets_data, planets_data, (Nplanets*Nparams)*sizeof(*planets->planets_data));
    memcpy(planets->masses, masses, (Nplanets+1)*sizeof(*planets->masses));

    return;
}

void rebx_free_planets_pointers(struct rebx_planets* const planets){
    free(planets->planets_data); 
    free(planets->masses); 

    return;
}

void rebx_free_planets(struct rebx_planets* const planets){
    rebx_free_planets_pointers(planets);
    free(planets);

    return;
}