/**********************************************************************
 *
 * Vector Space System Project / Material Library Module / IAPWS-IF97
 *
 * Copyright (C) 2020 CIAE.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 * File: property.v
 * Author: YU YANG
 * Created Time: 2020-05-07
 * Version: 0.0.1
 *
 **********************************************************************/
module if97
import math

//****************************************************************
//Calculate surface tension between water phase and vapour phase
//The correlation is based on the IAPWS Release on Surface Tension of Ordinary Water Substance, September 1994.
//****************************************************************
fn surface_tensions(T f64) f64{
	//temperature (K); surface tension (N/m)
	tau := 1.0 - T / if97_tcrit
	return 0.2358 * math.pow(tau, 1.256) * (1 - 0.625 * tau)
}

//****************************************************************
// Based on IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance
//****************************************************************
const (
	if97_mucrit = 1.0e-6		/* Pa-s */
    hk_data = [1.67752, 2.20462, 0.6366564, -0.241605]
	/*hij_data = [	ijn(0, 0, 5.20094e-1),
					ijn(1, 0, 8.50895e-2),
					ijn(2, 0, -1.08374),
					ijn(3, 0, -2.89555e-1),
					ijn(0, 1, 2.22531e-1),
					ijn(1, 1, 9.99115e-1),
					ijn(2, 1, 1.88797),
					ijn(3, 1, 1.26613),
					ijn(5, 1, 1.20573e-1),
					ijn(0, 2, -2.81378e-1),
					ijn(1, 2, -9.06851e-1),
					ijn(2, 2, -7.72479e-1),
					ijn(3, 2, -4.89837e-1),
					ijn(4, 2, -2.57040e-1),
					ijn(0, 3, 1.61913e-1),
					ijn(1, 3, 2.57399e-1),
					ijn(0, 4, -3.25372e-2),
					ijn(3, 4, 6.98452e-2),
					ijn(4, 5, 8.72102e-3),
					ijn(3, 6, -4.35673e-3),
					ijn(5, 6, -5.93264e-4)]
	hij_data=[ [5.20094E-1, 2.22531E-1, -2.81378E-1, 1.61913E-1, -3.25372E-2], 
				[8.50895E-2, 9.99115E-1, -9.06851E-1, 2.57399E-1], 
				[-1.08374, 1.88797, -7.72479E-1], 
				[-2.89555E-1, 1.26613, -4.89837E-1, 0.0, 6.98452E-2, 0.0, -4.35673E-3], 
				[0.0, 0.0, -2.57040E-1, 0.0, 0.0, 8.72102E-3], 
				[0.0, 1.20573E-1, 0.0, 0.0, 0.0, 0.0, -5.93264E-4]]
	*/
	hij_data=[ [5.20094E-1, 2.22531E-1, -2.81378E-1, 1.61913E-1, -3.25372E-2], 
				[8.50895E-2, 9.99115E-1, -9.06851E-1, 2.57399E-1], 
				[-1.08374, 1.88797, -7.72479E-1], 
				[-2.89555E-1, 1.26613, -4.89837E-1, 0.0, 6.98452E-2, 0.0, -4.35673E-3], 
				[-2.57040E-1, 8.72102E-3], 
				[1.20573E-1, -5.93264E-4]]
)

fn mu0(tau f64) f64{
    // viscosity in the dilute-gas limit
    mut sum := 0.0
    for i := 0; i < 4; i++ {
        sum += hk_data[i] * ipow(tau, i)
    }
    return 100.0 / (math.sqrt(tau) * sum)
}

/*	ijn() version
fn mu1(del f64, tau f64) f64{
    // contribution to viscosity due to finite density
    mut sum := 0.0
	for d in hij_data{
		sum += d.n_ * ipow(tau - 1.0, d.i_) * ipow(del - 1.0, d.j_)
	}
    return math.exp(del * sum)
}*/

// fast version of mu1()
fn mu1_1(del f64, tau f64) f64{
    // contribution to viscosity due to finite density
	para_i := tau - 1.0
	para_j := del -1.0
	mut sumi := 0.0
	mut sumj := 0.0
	for i := 5; i >= 0; i--{
		sumi *= para_i
		sumj = hij_data[i][hij_data[i].len-1]
		if i < 4 {
			for j := hij_data[i].len-2; j >= 0; j-- {
				sumj *= para_j
				sumj += hij_data[i][j]
			}
		}
		else {
			sumj *= ipow(para_j, 2*i-5)
			sumj += hij_data[i][0]
			sumj *= ipow(para_j, 6 - i)
		}
		sumi += sumj
	}
    return math.exp(del * sumi)
}

/* // another version of mu1()
fn mu1_2(del f64, tau f64) f64{
    // contribution to viscosity due to finite density
	mut sumi := 0.0
	mut sumj := 0.0
	para_i := tau - 1.0
	para_j := del -1.0
	for i := 5; i >= 0; i--{
		sumi *= para_i
		sumj = hij_data[i][hij_data[i].len-1]
		for j := hij_data[i].len-2; j >= 0; j-- {
			sumj *= para_j
			sumj += hij_data[i][j]
		}
		sumi += sumj
	}
    return math.exp(del * sumi)
}*/

fn mu2(del f64, tau f64) f64{
	// critical enhancement to viscosity
    return 1.0
}

pub fn mu_r1208(rho f64, T f64) f64{
	del := rho / if97_rhocrit
	tau := if97_tcrit / T

	mu2 := 1.0

	return if97_mucrit * mu0(tau) * mu1_1(del,tau) * mu2
}


//****************************************************************
// Thermal Conductivity [W/m.K],  Version1
// IAPWS Formulation 2011 for the Thermal Conductivity of Ordinary Water Substance
//****************************************************************
const (
	lk_data = [	2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4]
	lij_data = [	ijn(0, 0, 1.60397357),
					ijn(0, 1, -6.46013523e-1),
					ijn(0, 2, 1.11443906e-1),
					ijn(0, 3, 1.02997357e-1),
					ijn(0, 4, -5.04123634e-2),
					ijn(0, 5, 6.09859258e-3),
					ijn(1, 0, 2.33771842),
					ijn(1, 1, -2.78843778),
					ijn(1, 2, 1.53616167),
					ijn(1, 3, -4.63045512e-1),
					ijn(1, 4, 8.32827019e-2),
					ijn(1, 5, -7.19201245e-3),
					ijn(2, 0, 2.19650529),
					ijn(2, 1, -4.54580785),
					ijn(2, 2, 3.55777244),
					ijn(2, 3, -1.40944978),
					ijn(2, 4, 2.75418278e-1),
					ijn(2, 5, -2.05938816e-2),
					ijn(3, 0, -1.21051378),
					ijn(3, 1, 1.60812989),
					ijn(3, 2, -6.21178141e-1),
					ijn(3, 3, 7.16373224e-2),
					ijn(4, 0, -2.7203370),
					ijn(4, 1, 4.57586331),
					ijn(4, 2, -3.18369245),
					ijn(4, 3, 1.1168348),
					ijn(4, 4, -1.9268305e-1),
					ijn(4, 5, 1.2913842e-2)]
)

fn lamda0(tau f64) f64{
	//thermal conductivity in the dilute-gas limit
	mut sum := 0.0
    for i := 0; i < 5; i++ {
        sum += lk_data[i] * ipow(tau, i)
    }
    return 1.0 / (math.sqrt(tau) * sum)
}

fn lamda1(del f64, tau f64) f64{
	//contribution to thermal conductivity due to finite density
    mut sum := 0.0
	for d in lij_data{
		sum += d.n_ * ipow(tau - 1.0, d.i_) * ipow(del - 1, d.j_)
	}
    return math.exp(del * sum)
}

fn lamda2(del f64, tau f64) f64{
	return 0.0
}

pub fn k_r1511(rho f64, T f64) f64{
	del := rho / if97_rhocrit
	tau := if97_tcrit / T

	return lamda0(tau) * lamda1(del, tau) + lamda2(del, tau)
} 

// Thermal Conductivity [W/m.K],  Version2
pub fn k_rhot(rho f64, T f64) f64{
	tstar := 647.26
	rhostar := 317.7
	kstar := 1.0

	b0 := -0.397070
	b1 := 0.400302
	b2 := 1.060000
	bb1 := -0.171587
	bb2 := 2.392190

	c1 := 0.642857
	c2 := -4.11717
	c3 := -6.17937
	c4 := 0.00308976
	c5 := 0.0822994
	c6 := 10.0932

	d1 := 0.0701309
	d2 := 0.0118520
	d3 := 0.00169937
	d4 := -1.0200

	a_count := 4
	a := [0.0102811, 0.0299621, 0.0156146, -0.00422464]

	t_bar := T / tstar
	rhobar := rho / rhostar

	mut t_root := math.sqrt(t_bar)

	mut t_pow := t_root
	mut lam := 0.0
	for k := 0; k < a_count; k++ {
		lam += a[k] * t_pow
		t_pow *= t_bar
	}

	lam += b0 + b1 * rhobar + b2 * math.exp(bb1 * sq(rhobar + bb2))

	dt_bar := fabs(t_bar - 1) + c4
	dt_barpow := math.pow(dt_bar, 3.0/5)
	q := 2.0 + c5 / dt_barpow

	mut s := 0.0
	if t_bar >= 1 {
		s = 1.0 / dt_bar
	}else{
		s = c6 / dt_barpow
	}

	rhobar18 := math.pow(rhobar, 1.8)
	rhobar_q := math.pow(rhobar, q)

	lam += (d1 / ipow(t_bar,10) + d2) * rhobar18 * math.exp(c1 * (1 - rhobar * rhobar18)) + d3 * s * rhobar_q *	math.exp((q/(1+q))*(1 - rhobar*rhobar_q)) + d4 * math.exp(c2 * ipow(t_root,3) + c3 / ipow(rhobar,5))

	return kstar * lam
}
