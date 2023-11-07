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

int rc_idx(const int r, const int c, const int N);

void rebx_planets_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    const double t = sim->t;
    const double G = sim->G;

    // planets data
    struct rebx_planets* const planets = rebx_get_param(rebx, force->ap, "planets_data");
    const int Nplanets = planets->Nplanets;
    const int Nparams = planets->Nparams;
    const double* planets_data = planets->planets_data;
    const double* masses = planets->masses;
    const double star_mass = masses[0];

    // printf("\n");
    // for (int i=0; i<Nplanets; i++){
    //     for (int j=0; j<Nparams; j++){
    //         printf("%e ", planets_data[rc_idx(i, j, Nparams)]);
    //     }
    // }

    // calculate star and planet positions at time t
    double xyz[Nplanets + 1][3];
    xyz[0][0] = 0.;
    xyz[0][1] = 0.;
    xyz[0][2] = 0.;

    for (int i=0; i<Nplanets; i++) {
        double a     = planets_data[rc_idx(i, 0, Nparams)];
        double e     = planets_data[rc_idx(i, 1, Nparams)];
        double inc   = planets_data[rc_idx(i, 2, Nparams)];
        double omega = planets_data[rc_idx(i, 3, Nparams)];
        double Omega = planets_data[rc_idx(i, 4, Nparams)];
        double M0    = planets_data[rc_idx(i, 5, Nparams)];
        double n     = planets_data[rc_idx(i, 6, Nparams)];

        // calculate the mean motion at time t
        double M = M0 + n * t;
        // convert to true anomaly
        double f = reb_tools_M_to_f(e, M);

        // find distance from barycenter at true anomaly f
        double rp = a * (1. - e*e) / (1. + e*cos(f));
        // get xyz position of planet
        xyz[i+1][0] = rp*(cos(Omega)*cos(omega+f) - sin(Omega)*sin(omega+f)*cos(inc));
        xyz[i+1][1] = rp*(sin(Omega)*cos(omega+f) + cos(Omega)*sin(omega+f)*cos(inc));
        xyz[i+1][2] = rp*(sin(inc)*sin(omega+f));

        // star's COM position
        xyz[0][0] += -masses[i+1]*xyz[i+1][0]/star_mass;
        xyz[0][1] += -masses[i+1]*xyz[i+1][1]/star_mass;
        xyz[0][2] += -masses[i+1]*xyz[i+1][2]/star_mass;
    }
    // printf("%e %e %e\n", xyz[0][0], xyz[0][1], xyz[0][2]);

    // calculate force from planet/star to particles
    for (int i=0; i<Nplanets + 1; i++){
        double GM = -(G * masses[i]);
        double x = xyz[i][0];
        double y = xyz[i][1];
        double z = xyz[i][2];

        for (int j=0; j<N; j++){
            // get distance to test particle from massive particle
            double dx = particles[j].x - x;
            double dy = particles[j].y - y;
            double dz = particles[j].z - z;
            double d1 = sqrt(dx*dx + dy*dy + dz*dz);
            double d3 = d1*d1*d1;

            // add acceleration from massive particle to test particle
            particles[j].ax += GM * dx/d3 ;
            particles[j].ay += GM * dy/d3 ;
            particles[j].az += GM * dz/d3 ;
        }
    }

    
}

int rc_idx(const int r, const int c, const int N){
    // index of a two dimensional array, given row r and column c
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