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
 * pf_ap (double)               Yes         Semimajor axis of the planet.
 * pf_mplanet (double)          Yes         Planet mass.
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

    // get particles and constants
    const double* ap = rebx_get_param(rebx, force->ap, "pf_ap");
    const double* inc = rebx_get_param(rebx, force->ap, "pf_inc");
    const double* m_planet = rebx_get_param(rebx, force->ap, "pf_mplanet");
    const double mu = sim->G * (particles[0].m + *m_planet);
    const double t = sim->t;

    // mean motion * t
    double npt = sqrt( mu / pow(*ap, 3) ) * t;

    // get planet position
    double planet_x = (*ap) * cos(npt);
    double planet_y = (*ap) * sin(npt) * cos(*inc);
    double planet_z = (*ap) * sin(npt) * sin(*inc);

    for (int i=0; i<N; i++){
        if (i > 0){
            // get particle position relative to planet
            double dx = (particles[i].x - planet_x);
            double dy = (particles[i].y - planet_y);
            double dz = (particles[i].z - planet_z);
            double dist = sqrt( pow(dx, 2) + pow(dy, 2) + pow(dz, 2) );
            double d3 = pow(dist, 3);

            // update force
            particles[i].ax += *m_planet/d3 * dx;
            particles[i].ay += *m_planet/d3 * dy;
            particles[i].az += *m_planet/d3 * dz;
        }
    }
}
