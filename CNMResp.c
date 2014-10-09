void CNMResp (estate_type *es, cstate_type *cs, pstate_type *ps, cflux_type *cf, epc_type *con, int nlevgrnd)
{
/*
 * DESCRIPTION:
 * Module holding maintenance respiration routines for coupled carbon
 * nitrogen code.
 */

    int c,p,j;          /* indices */
    double br;          /* base rate (gC/gN/s) */
    double q10;         /* temperature dependence */
    double tc;          /* temperature correction, 2m air temp (unitless) */
    double tcsoi[nlevgrnd];      /* temperature correction by soil layer (unitless) */

    /* base rate for maintenance respiration is from:
     * M. Ryan, 1991. Effects of climate change on plant respiration.
     * Ecological Applications, 1(2), 157-167.
     * Original expression is br = 0.0106 molC/(molN h)
     * Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)*/

    br = 2.525e-6;

    /* Peter Thornton: 3/13/09
     * Q10 was originally set to 2.0, an arbitrary choice,
     * but reduced to 1.5 as part of the tuning
     * to improve seasonal cycle of atmospheric CO2 concentration in global
     * simulatoins */

    q10 = 1.5;

    /* column loop to calculate temperature factors in each soil layer */
    for (j = 0; j < nlevgrnd; j++)
    {
        /* calculate temperature corrections for each soil layer, for use in
         * estimating fine root maintenance respiration with depth */

         tcsoi[j] = pow(q10, (es->t_soisno[j] - SHR_CONST_TKFRZ - 20.0) / 10.0);
    }

    /* pft loop for leaves and live wood */
    /* calculate maintenance respiration fluxes in
     * gC/m2/s for each of the live plant tissues.
     * Leaf and live wood MR */

    tc = pow(q10, (es->t_ref2m - SHR_CONST_TKFRZ - 20.0) / 10.0);
    cf->leaf_mr = ns->leafn * br * tc;
    if (con->woody == 1)
    {
        cf->livestem_mr = ns->livestemn * br * tc;
        cf->livecroot_mr = ns->livecrootn * br * tc;
    }

    /* soil and pft loop for fine root */
    for (j = 0; j < nlevgrnd; j++)
    {
        /* Fine root MR
         * rootfr(j) sums to 1.0 over all soil layers, and
         * describes the fraction of root mass that is in each
         * layer.  This is used with the layer temperature correction
         * to estimate the total fine root maintenance respiration as a
         * function of temperature and N content. */

         cf->froot_mr = cf->froot_mr + ns->frootn * br * tcsoi[j] * ps->rootfr[j];
    }
}
