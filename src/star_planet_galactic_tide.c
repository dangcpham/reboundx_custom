/** * @file star_planet_galactic_tide.c
 * @brief   Forces from star, planet, and galactic tide
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
 * spgt_inc (double)              Yes         Inclination of the planet.
 * spgt_ap (double)               Yes         Semi-major axis of the planet.
 * spgt_as (double)               Yes         Semi-major axis of the star.
 * spgt_n (double)                Yes          Mean motion of the planet.
 * spgt_m0p (double)              Yes          Initial planet mean anomaly.
 * spgt_mplanet (double)          Yes         Planet mass.
 * spgt_mstar (double)            Yes         Star mass.
 * spgt_rho (double)              Yes         Galactic density.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_spgt(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;

    // get physical constants and particles
    const double t = sim->t;
    const double G = sim->G;

    // get constants
    const double* inc_p = rebx_get_param(rebx, force->ap, "spgt_inc");
    const double* ap = rebx_get_param(rebx, force->ap, "spgt_ap");
    const double* as = rebx_get_param(rebx, force->ap, "spgt_as");
    const double* n = rebx_get_param(rebx, force->ap, "spgt_n");
    const double* m0p = rebx_get_param(rebx, force->ap, "spgt_m0p");
    const double* m_planet = rebx_get_param(rebx, force->ap, "spgt_mplanet");
    const double* m_star = rebx_get_param(rebx, force->ap, "spgt_mstar");
    const double* rho = rebx_get_param(rebx, particles[0].ap, "spgt_rho");

    // other constants
    double inc_s = -1. * (*inc_p);
    double nt = (*n) * t + (*m0p);
    double Gm_p = G * (*m_planet);
    double Gm_s = G * (*m_star);
    double Gm_com = G * (*m_planet + *m_star);

    // particle position
    double x = particles[0].x;
    double y = particles[0].y;
    double z = particles[0].z;

    // distance from COM
    double dcom = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );    

    // galactic tide
    double azgt = -4.0 * M_PI * (*rho) * G * z;

    if (dcom > 1000.) {
        double d3com = pow(dcom,3);
        double axcom = (-Gm_com / d3com) * x;
        double aycom = (-Gm_com / d3com) * y;
        double azcom = (-Gm_com / d3com) * z;        

        particles[0].ax += axcom;
        particles[0].ay += aycom;
        particles[0].az += azcom + azgt;
    } else {
        // get planet position
        double planet_x = (*ap) * cos(nt);
        double planet_y = (*ap) * sin(nt) * cos(*inc_p);
        double planet_z = (*ap) * sin(nt) * sin(*inc_p);

        // get star position
        double star_x = -1. * (*as) * cos(nt);
        double star_y = -1. * (*as) * sin(nt) * cos(inc_s);
        double star_z = (*as) * sin(nt) * sin(inc_s);

        // get comet distance to planet
        double dxp = (x - planet_x);
        double dyp = (y - planet_y);
        double dzp = (z - planet_z);
        double dp = sqrt( pow(dxp,2) + pow(dyp,2) + pow(dzp,2) );
        double d3p = pow(dp,3);

        // get star distance to planet
        double dxs = (x - star_x);
        double dys = (y - star_y);
        double dzs = (z - star_z);
        double ds = sqrt( pow(dxs,2) + pow(dys,2) + pow(dzs,2) );
        double d3s = pow(ds,3);

        // calculate force
        double axp = (-Gm_p/d3p) * dxp;
        double ayp = (-Gm_p/d3p) * dyp;
        double azp = (-Gm_p/d3p) * dzp;

        double axs = (-Gm_s/d3s) * dxs;
        double ays = (-Gm_s/d3s) * dys;
        double azs = (-Gm_s/d3s) * dzs;

        // update force
        particles[0].ax += axp + axs;
        particles[0].ay += ayp + ays;
        particles[0].az += azp + azs + azgt;
    }
}