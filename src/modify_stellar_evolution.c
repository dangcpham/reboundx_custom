#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_modify_stellar_evolution(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    struct rebx_interpolator* stellarinterp = rebx_get_param(sim->extras, sim->particles[0].ap, "stellar_interp");
    sim->particles[0].m = rebx_interpolate(sim->extras, stellarinterp, sim->t);
    reb_simulation_move_to_com(sim);
    // printf("Stellar mass: %.17g\n", sim->particles[0].m);
}