/** * @file planet_force.c
 * @brief   Force due to flybys for massless particles
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
 * pf_inc (double)              Yes         Inclination of the planet.
 * pf_ap (double)               Yes         Semi-major axis of the planet.
 * pf_as (double)               Yes         Semi-major axis of the star.
 * pf_n (double)               Yes          Mean motion of the planet.
 * pf_m0p (double)             Yes          Initial planet mean anomaly.
 * pf_mplanet (double)          Yes         Planet mass.
 * pf_mstar (double)            Yes         Star mass.
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

// void rebx_xyz_unbound(struct rebx_extras* const rebx, const double* t0, const double* tf, const double* a, const double* ecc, const double* E0, const double* Omega, const double* inc, const double* omega);

struct rebx_flybys* rebx_create_flybys(struct rebx_extras* const rebx, const int Nvalues, const double* t0, const double* tf, const double* a, const double* ecc, const double* E0, const double* Omega, const double* inc, const double* omega){
    struct rebx_flybys* flybys = rebx_malloc(rebx, sizeof(*flybys));
    rebx_init_flybys(rebx, flybys, Nvalues, t0, tf, a, ecc, E0, Omega, inc, omega);
    return flybys;
}

void rebx_init_flybys(struct rebx_extras* const rebx, struct rebx_flybys* const flybys, const int Nvalues, const double* t0, const double* tf, const double* a, const double* ecc, const double* E0, const double* Omega, const double* inc, const double* omega){
    flybys->Nvalues = Nvalues;
    flybys->t0 = calloc(Nvalues, sizeof(*flybys->t0));
    flybys->tf = calloc(Nvalues, sizeof(*flybys->tf));
    flybys->a  = calloc(Nvalues, sizeof(*flybys->a));
    flybys->ecc = calloc(Nvalues, sizeof(*flybys->ecc));
    flybys->E0 = calloc(Nvalues, sizeof(*flybys->E0));
    flybys->Omega = calloc(Nvalues, sizeof(*flybys->Omega));
    flybys->inc = calloc(Nvalues, sizeof(*flybys->inc));
    flybys->omega = calloc(Nvalues, sizeof(*flybys->omega));

    memcpy(flybys->t0, t0, Nvalues*sizeof(*flybys->t0));
    memcpy(flybys->tf, tf, Nvalues*sizeof(*flybys->tf));
    memcpy(flybys->a, a, Nvalues*sizeof(*flybys->a));
    memcpy(flybys->ecc, ecc,  Nvalues*sizeof(*flybys->ecc));
    memcpy(flybys->E0, E0,  Nvalues*sizeof(*flybys->E0));
    memcpy(flybys->Omega, Omega, Nvalues*sizeof(*flybys->Omega));
    memcpy(flybys->inc, inc, Nvalues*sizeof(*flybys->inc));
    memcpy(flybys->omega, omega, Nvalues*sizeof(*flybys->omega));

    return;
}

void rebx_free_flybys_pointers(struct rebx_flybys* const flybys){
    free(flybys->t0); 
    free(flybys->tf);
    free(flybys->a);
    free(flybys->ecc);
    free(flybys->E0);
    free(flybys->Omega);
    free(flybys->inc);
    free(flybys->omega);

    return;
}

void rebx_free_flybys(struct rebx_flybys* const flybys){
    rebx_free_flybys_pointers(flybys);
    free(flybys);

    return;
}