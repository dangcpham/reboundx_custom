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

void rebx_xyz_unbound(double* x, double* y, double* z, double a, double ecc, double ecc_anomaly, double Omega, double inc, double omega){
    // calculates position of flyby star relative to the COM of the current sim
    // e.g. if there's only star in the sim, then xyz is the heliocentric coordinates of the flyby star

    double cosh_term = (ecc-cosh(ecc_anomaly));
    double sqrte_term = sqrt(ecc*ecc - 1);

    *x = -a*(sinh(ecc_anomaly)*sqrte_term*(sin(omega)*cos(Omega) + sin(Omega)*cos(inc)*cos(omega)) -
             cosh_term*(sin(Omega)*sin(omega)*cos(inc) - cos(Omega)*cos(omega)));

    *y = -a*(sinh(ecc_anomaly)*sqrte_term*(sin(Omega)*sin(omega) - cos(inc)*cos(Omega)*cos(omega)) +
             cosh_term*(sin(omega)*cos(inc)*cos(Omega) + sin(Omega)*cos(omega)));

    *z = -a*sin(inc)*(cosh_term*sin(omega)-sinh(ecc_anomaly)*cos(omega)*sqrte_term);

}

void rebx_flybys_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    struct rebx_flybys* const flybys = rebx_get_param(rebx, force->ap, "flybys_data");
    const double Nflybys = flybys->Nvalues;
    const double t = sim->t;
    const double G = sim->G;

    for (int i=0; i<Nflybys; i++){
        if ((t >= flybys->t0[i]) && (t <= flybys->tf[i])){
            // calculate new eccentric anomaly at time t
            const double M = (2.*M_PI - flybys->M0[i] + (flybys->n[i] * t));
            const double ecc_anomaly = reb_tools_M_to_E(flybys->ecc[i], M);

            // calculate position of flyby at time t
            double x, y, z;
            rebx_xyz_unbound(&x, &y, &z, flybys->a[i], flybys->ecc[i], ecc_anomaly, flybys->Omega[i], flybys->inc[i], flybys->omega[i]);

            // calculate force from x, y, z to each particle
            double GM = G * flybys->m[i];
            for (int j=0; j<N; j++){
                // ignore particles that don't feel the flyby force
                // e.g. this would be the central star if there is a star + massless test particle
                // because flyby implemented is in heliocentric coordinate system
                const int* flybys_mode = rebx_get_param(rebx, particles[j].ap, "flybys_mode");

                // printf("%f %i %i %f %f %f %f\n ", flybys->ecc[i], j, *flybys_mode, t, x, y, z);
                if (*flybys_mode == 1.0){
                    // get particle distance to flyby star
                    double dx = (particles[j].x - x);
                    double dy = (particles[j].y - y);
                    double dz = (particles[j].z - z);
                    double d = sqrt( dx*dx + dy*dy + dz*dz );
                    double d3 = d*d*d;

                    // update particle acceleration
                    particles[j].ax -= (GM/d3) * dx ;
                    particles[j].ay -= (GM/d3) * dy ;
                    particles[j].az -= (GM/d3) * dz ;
                }
            }
        }
    }
}

struct rebx_flybys* rebx_create_flybys(struct rebx_extras* const rebx, const int Nvalues, const double* m, const double* t0, const double* tf, const double* n, const double* a, const double* ecc, const double* M0, const double* Omega, const double* inc, const double* omega){
    struct rebx_flybys* flybys = rebx_malloc(rebx, sizeof(*flybys));
    rebx_init_flybys(rebx, flybys, Nvalues, m, t0, tf, n, a, ecc, M0, Omega, inc, omega);
    return flybys;
}

void rebx_init_flybys(struct rebx_extras* const rebx, struct rebx_flybys* const flybys, const int Nvalues, const double* m, const double* t0, const double* tf, const double* n, const double* a, const double* ecc, const double* M0, const double* Omega, const double* inc, const double* omega){
    flybys->Nvalues = Nvalues;
    flybys->m = calloc(Nvalues, sizeof(*flybys->m));
    flybys->t0 = calloc(Nvalues, sizeof(*flybys->t0));
    flybys->tf = calloc(Nvalues, sizeof(*flybys->tf));
    flybys->n = calloc(Nvalues, sizeof(*flybys->n));
    flybys->a = calloc(Nvalues, sizeof(*flybys->a));
    flybys->ecc = calloc(Nvalues, sizeof(*flybys->ecc));
    flybys->M0 = calloc(Nvalues, sizeof(*flybys->M0));
    flybys->Omega = calloc(Nvalues, sizeof(*flybys->Omega));
    flybys->inc = calloc(Nvalues, sizeof(*flybys->inc));
    flybys->omega = calloc(Nvalues, sizeof(*flybys->omega));

    memcpy(flybys->m, m, Nvalues*sizeof(*flybys->m));
    memcpy(flybys->t0, t0, Nvalues*sizeof(*flybys->t0));
    memcpy(flybys->tf, tf, Nvalues*sizeof(*flybys->tf));
    memcpy(flybys->n, n, Nvalues*sizeof(*flybys->n));
    memcpy(flybys->a, a, Nvalues*sizeof(*flybys->a));
    memcpy(flybys->ecc, ecc,  Nvalues*sizeof(*flybys->ecc));
    memcpy(flybys->M0, M0,  Nvalues*sizeof(*flybys->M0));
    memcpy(flybys->Omega, Omega, Nvalues*sizeof(*flybys->Omega));
    memcpy(flybys->inc, inc, Nvalues*sizeof(*flybys->inc));
    memcpy(flybys->omega, omega, Nvalues*sizeof(*flybys->omega));

    return;
}

void rebx_free_flybys_pointers(struct rebx_flybys* const flybys){
    free(flybys->m); 
    free(flybys->t0); 
    free(flybys->tf);
    free(flybys->n);
    free(flybys->a);
    free(flybys->ecc);
    free(flybys->M0);
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