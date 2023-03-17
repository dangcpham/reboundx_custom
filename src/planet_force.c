/** * @file planet_force.c
 * @brief   Force due to planet for massless particles
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
 * pf_mplanet (double)          Yes         Planet mass.
 * pf_mstar (double)            Yes         Star mass.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_planet_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;

    // get physical constants and particles
    const double t = sim->t;
    const double G = sim->G;

    // get particle setups
    const double* inc_p = rebx_get_param(rebx, force->ap, "pf_inc");
    const double* ap = rebx_get_param(rebx, force->ap, "pf_ap");
    const double* as = rebx_get_param(rebx, force->ap, "pf_as");
    const double* n = rebx_get_param(rebx, force->ap, "pf_n");
    const double* m_planet = rebx_get_param(rebx, force->ap, "pf_mplanet");
    const double* m_star = rebx_get_param(rebx, force->ap, "pf_mstar");

    // other constants
    double inc_s = -1. * (*inc_p);
    double nt = (*n) * t;
    double Gm_p = G * (*m_planet);
    double Gm_s = G * (*m_star);

    // particle position
    double x = particles[0].x;
    double y = particles[0].y;
    double z = particles[0].z;

    for (int i=0; i<N; i++){
        // force from planet
        double ax_p, ay_p, az_p;
        force_from_planet(nt, *inc_p, *ap, Gm_p, x, y, z, &ax_p, &ay_p, &az_p);

        // force from star
        double ax_s, ay_s, az_s;
        force_from_star(  nt,  inc_s, *as, Gm_s, x, y, z, &ax_s, &ay_s, &az_s);

        // update force
        particles[i].ax += ax_p + ax_s;
        particles[i].ay += ay_p + ay_s;
        particles[i].az += az_p + az_s;
    }
}

void force_from_planet(const double npt, const double inc, const double a, const double Gm, const double px, const double py, const double pz, double* ax, double* ay, double* az){
    // get massive particle position
    double planet_x = a * cos(npt);
    double planet_y = a * sin(npt) * cos(inc);
    double planet_z = a * sin(npt) * sin(inc);

    // get massless position relative to massive
    double dx = (px - planet_x);
    double dy = (py - planet_y);
    double dz = (pz - planet_z);
    double dist = sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2) );
    double d3 = pow(dist,3);

    // update force from massive
    *ax = -Gm/d3 * dx;
    *ay = -Gm/d3 * dy;
    *az = -Gm/d3 * dz;
}

void force_from_star(const double npt, const double inc, const double a, const double Gm, const double px, const double py, const double pz, double* ax, double* ay, double* az){
    // get massive particle position
    double planet_x = -1. * a * cos(npt);
    double planet_y = -1. * a * sin(npt) * cos(inc);
    double planet_z = a * sin(npt) * sin(inc);

    // get massless position relative to massive
    double dx = (px - planet_x);
    double dy = (py - planet_y);
    double dz = (pz - planet_z);
    double dist = sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2) );
    double d3 = pow(dist,3);

    // update force from massive
    *ax = -Gm/d3 * dx;
    *ay = -Gm/d3 * dy;
    *az = -Gm/d3 * dz;
}