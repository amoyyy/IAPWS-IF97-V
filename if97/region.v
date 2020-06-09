/**********************************************************************
 *
 * Vector Space System Project / Material Library Module / IF97
 *
 * Copyright (C) 2020 CIAE.
 *
 * This is free software you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 * File: region.v
 * Author: YU YANG
 * Created Time: 2020-05-07
 * Version: 0.0.1
 *
 **********************************************************************/
module if97
import math

/*
    * v: Specific volume, [m³/kg]
    * h: Specific enthalpy, [kJ/kg]
    * s: Specific entropy, [kJ/kgK]
    * cp: Specific isobaric heat capacity, [kJ/kgK]
    * cv: Specific isocoric heat capacity, [kJ/kgK]
    * w: Speed of sound, [m/s]
    * alphav: Cubic expansion coefficient, [1/K]
    * kt: Isothermal compressibility, [1/MPa]
*/

//****************************************************************
// REGION 1 G(p,T) EQUATIONS, IF97-Rev, Eq 7,
// 273.15 K < T < 623.15 K, ps(T) < p < 100 MPa.
//****************************************************************
// Internal function
fn r1_pitau(p f64, t f64) (f64, f64) {
	return 7.1 - p / pstar_region1, tstar_region1/ t -1.222
}

// Public function
pub fn r1_v(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	k1r := (7.1 - pii) / pii
	mut gamp := 0.0
	for d in pt_r1 {
		gamp += -d.n_ * d.i_ * ipow2(pii, taui, d.i_, d.j_)
	}
	return if97_r *  T  * gamp * k1r / p / 1000.0
}

pub fn r1_u(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	k1 := pii / (7.1 - pii)
	k2 := taui / (taui + 1.222) 
	mut sum := 0.0
	for d in pt_r1 {
		coeff := d.j_ / k2 + d.i_ / k1
		sum += coeff * d.n_ * ipow2(pii, taui, d.i_, d.j_)
	}
	return if97_r *  T * sum
}

pub fn r1_s(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	k2 := taui / (taui + 1.222) 
	mut sum := 0.0
	for d in pt_r1 {
		coeff := d.j_ / k2 - 1.0
		sum += coeff * d.n_ * ipow2(pii, taui, d.i_, d.j_)
	}
	return if97_r *  sum
}

pub fn r1_h(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	k2r := (taui + 1.222) / taui
	mut gamt := 0.0
	for d in pt_r1 {
		gamt += d.n_ * d.j_ * ipow2(pii, taui, d.i_, d.j_)
	}
	return if97_r *  T * gamt * k2r 
}

pub fn r1_cp(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	tau := (taui + 1.222)
	k2rs := sq(tau/taui)
	mut gamtt := 0.0
	for d in pt_r1 {
		gamtt += d.n_ * d.j_ * (d.j_ - 1) * ipow2(pii, taui, d.i_, d.j_)
	}
	return -(if97_r *  k2rs * gamtt)
}

pub fn r1_cv(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	tau := (taui + 1.222)
	mut gamtt := 0.0
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in pt_r1 {
		tmp := d.n_ * ipow2(pii, taui, d.i_-2, d.j_-2)
		gamtt += d.j_ * (d.j_ - 1) * tmp
		sum1 += d.i_ * (tau * d.j_ - taui) * tmp
		sum2 += d.i_ * (d.i_ - 1) * tmp
	}
	gamtt *= sq(pii)
	sum1 *= pii
	return if97_r *  (-sq(tau)*gamtt + sq(sum1)/sum2)
}

pub fn r1_w(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	k0 := pii * taui
	k2 := taui / (taui + 1.222) 
	mut sum1 := 0.0
	mut sum2 := 0.0
	mut sum3 := 0.0
	mut sum4 := 0.0
	for d in pt_r1 {
		tmp := d.n_ * ipow2(pii, taui, d.i_-1, d.j_-1)
		tmp1 := d.i_ * tmp
		sum1 += tmp1
		sum2 += (d.j_ - k2) * tmp1
		sum3 += (d.j_ - 1) * d.j_ * tmp
		sum4 += (d.i_ - 1) * tmp1
	}
	return fabs(sum1) * math.sqrt(1000.0 * if97_r *  T * k0 / (sq(sum2)/sum3 - sum4))
}

pub fn r1_g(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	mut gam := 0.0
	for d in pt_r1 {
		gam += d.n_ * ipow2(pii, taui, d.i_, d.j_)
	}
	return if97_r *  T * gam
}

pub fn r1_a(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	pi := 7.1 - pii
	k1 := pi/pii
	mut sum := 0.0
	for d in pt_r1 {
		sum += d.n_ * (1.0 + d.i_* k1) * ipow2(pii, taui, d.i_, d.j_)
	}
	return if97_r *  T * sum
}

pub fn r1_alphav(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	k2r := (taui + 1.222) / taui
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in pt_r1 {
		tmp := d.n_ * d.i_ * ipow2(pii, taui, d.i_-1, d.j_-1)
		sum2 += tmp
		sum1 += tmp * d.j_
	}
	return (1.0 - k2r * sum1 / sum2) / T
}

pub fn r1_kt(p f64, T f64) f64{
	pii, taui := r1_pitau(p, T)
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in pt_r1 {
		tmp := d.n_ * d.i_ * ipow(pii, d.i_ - 2) * ipow(taui, d.j_)
		sum2 += tmp
		sum1 += tmp * (d.i_ - 1)
	}
	return sum1/sum2/pii/pstar_region1
}

//****************************************************************
// REGION 2 G(p,T) EQUATIONS, IF97-Rev, Eq 15-17,
// 273.15 K < T < 623.15 K, 0 < p < p_s(T).
// 623.15 K < T < 863.15 K, 0 < p < b23_p(T).
// 863.15 K < T < 1073.15 K, 0 < p < 100 MPa.
//****************************************************************
// Internal function
fn r2_pitau(p f64, t f64) (f64, f64) {
	return p / pstar_region245, tstar_region2/ t
}

// Public function
pub fn r2_v(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	mut gamp := 0.0
	for d in pt1_r2{
		gamp += d.n_ * d.i_ * ipow2(pi, taui, d.i_, d.j_)
	}
	gamp += 1.0
	return if97_r *  T  * gamp / p / 1000.0
}

pub fn r2_u(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	k2 := tau / taui
	mut sum := 0.0
	for d in pt0_r2{
		sum += d.n_ * d.j_ * ipow(tau, d.j_)
	}
	for d in pt1_r2{
		sum += d.n_ * (k2*d.j_ - d.i_) * ipow2(pi, taui, d.i_, d.j_) 
	}
	sum -= 1.0
	return if97_r *  T * sum
}

pub fn r2_s(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau - 0.5
	k2 := tau / taui
	mut sum := 0.0
	for d in pt0_r2{
		sum += d.n_ * (d.j_ - 1) * ipow(tau, d.j_)
	}
	for d in pt1_r2{
		sum += d.n_ * (k2*d.j_ - 1) * ipow2(pi, taui, d.i_, d.j_) 
	}
	sum -= math.log(pi)
	return if97_r *  sum
}

pub fn r2_h(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	k2 := tau / taui
	mut gamt0 := 0.0
	mut gamtr := 0.0
	for d in pt0_r2{
		gamt0 += d.n_ * d.j_ * ipow(tau, d.j_)
	}
	for d in pt1_r2{
		gamtr += d.n_ * d.j_ * ipow2(pi, taui, d.i_, d.j_)
	}
	return if97_r *  T * (gamt0 + k2 * gamtr)
}

pub fn r2_cp(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	k2s := sq(tau / taui)
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in pt0_r2{
		sum1 += d.n_ * d.j_ * (1 - d.j_) * ipow(tau, d.j_)
	}
	for d in pt1_r2{
		sum2 += d.n_ * d.j_ * (1 - d.j_) * ipow2(pi, taui, d.i_, d.j_) 
	}
	return if97_r *  (sum1 + k2s * sum2)
}

pub fn r2_cv(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	k2 := tau / taui
	k2s := sq(k2)
	
	mut ttg0tt := 0.0
	mut grtt := 0.0
	mut grpt := 0.0
	mut pgrp := 0.0
	mut ppgrpp := 0.0
	for d in pt0_r2{
		ttg0tt += d.n_ * d.j_* (d.j_ - 1) * ipow(tau, d.j_)
	}
	for d in pt1_r2{
		tmp := d.n_ * ipow2(pi, taui, d.i_, d.j_)
		pgrp += d.i_ * tmp
		grpt += d.i_ * d.j_ * tmp
		ppgrpp += d.i_ * (d.i_ - 1) * tmp
		grtt += d.j_ * (d.j_ - 1) * tmp
	}
	return if97_r *  (-(ttg0tt + k2s*grtt) - sq(1.0 + pgrp - k2*grpt) / (1.0 - ppgrpp))
}

pub fn r2_w(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	k2 := tau / taui
	k2s := sq(k2)

	mut ttg0tt := 0.0
	mut grtt := 0.0
	mut grpt := 0.0
	mut pgrp := 0.0
	mut ppgrpp := 0.0
	for d in pt0_r2{
		ttg0tt += d.n_ * d.j_* (d.j_ - 1) * ipow(tau, d.j_)
	}
	for d in pt1_r2{
		tmp := d.n_ * ipow2(pi, taui, d.i_, d.j_)
		pgrp += d.i_ * tmp
		grpt += d.i_ * d.j_ * tmp
		ppgrpp += d.i_ * (d.i_ - 1) * tmp
		grtt += d.j_ * (d.j_ - 1) * tmp
	}

	return math.sqrt(1000.0 * if97_r *  T * sq(1.0 + pgrp)/(1.0 - ppgrpp + sq(1.0 + pgrp - k2*grpt)/(ttg0tt + k2s*grtt)))
}

pub fn r2_g(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	mut gam := 0.0
	for d in pt0_r2{
		gam += d.n_ * ipow(tau, d.j_)
	}
	for d in pt1_r2{
		gam += d.n_ * ipow2(pi, taui, d.i_, d.j_)
	}
	gam += math.log(pi)
	return if97_r *  T * gam
}	

pub fn r2_a(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	mut sum := 0.0
	for d in pt0_r2{
		sum += d.n_ * ipow(tau, d.j_)
	}
	for d in pt1_r2{
		sum += d.n_ * (1 - d.i_) * ipow2(pi, taui, d.i_, d.j_)
	}
	sum += (math.log(pi) - 1.0)
	return if97_r *  T * sum
}

pub fn r2_alphav(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	k2 := tau / taui
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in pt1_r2{
		tmp := d.n_ * d.i_ * ipow2(pi, taui, d.i_, d.j_)
		sum2 += tmp
		sum1 += d.j_ * tmp
	}
	return (1.0 - k2 * sum1 / (1.0 + sum2)) / T
}

pub fn r2_kt(p f64, T f64) f64{
	pi, tau := r2_pitau(p, T)
	taui := tau-0.5
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in pt1_r2{
		tmp := d.n_ * d.i_ * ipow2(pi, taui, d.i_, d.j_)
		sum2 += tmp
		sum1 += (d.i_ - 1) * tmp
	}
	return (1.0-sum1) / (1.0+sum2) / p
}

//****************************************************************
// REGION 3 Helmholtz free energy f(rho,T) EQUATIONS, IF97-Rev, Eq 28
// 623.15 K < T < T(p){Eq 6} 
// p(T){Eq 5} < p < 100 MPa
//****************************************************************
// Internal function
fn r3_deltau(rho f64, t f64) (f64, f64) {
	return rho / rhostar_region3, tstar_region3 / t
}

// Public function
pub fn r3_p(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut dphid := 0.0
	for d in rt_r3 {
		dphid += d.n_ * d.i_ * ipow2(del, tau, d.i_, d.j_)
	}
	dphid += n1_r3 
	return rho * if97_r *  T * dphid / 1000.0
}

pub fn r3_u(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut tphit := 0.0
	for d in rt_r3 {
		tphit += d.n_ * d.j_ * ipow2(del, tau, d.i_, d.j_)
	}
	return if97_r  * tphit * T
}

pub fn r3_s(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut sum := 0.0
	for d in rt_r3 {
		sum += d.n_ * (d.j_ - 1)* ipow2(del, tau, d.i_, d.j_)
	}
	sum -= n1_r3*math.log(del)
	return if97_r *  sum
}

pub fn r3_h(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut sum := 0.0
	for d in rt_r3 {
		sum += d.n_ * (d.i_ + d.j_)* ipow2(del, tau, d.i_, d.j_)
	}
	sum += n1_r3
	return if97_r *  T * sum
}

pub fn r3_cp(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut ttphitt := 0.0
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in rt_r3 {
		tmp := d.n_ * ipow2(del, tau, d.i_, d.j_)
		ttphitt += d.j_ * (d.j_ - 1) * tmp
		sum1 += d.i_ * (1 - d.j_) * tmp
		sum2 += d.i_ * (d.i_ + 1) * tmp
	}
	sum1 += n1_r3
	sum2 += n1_r3
	return if97_r *  (sq(sum1) / sum2 - ttphitt)
}

pub fn r3_cv(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut ttphitt := 0.0
	for d in rt_r3 {
		ttphitt += d.n_ * d.j_ * (d.j_ - 1)* ipow2(del, tau, d.i_, d.j_)
	}
	return  - if97_r *  ttphitt 
}

pub fn r3_w(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut ttphitt := 0.0
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in rt_r3 {
		tmp := d.n_ * ipow2(del, tau, d.i_, d.j_)
		ttphitt += d.j_ * (d.j_ - 1) * tmp
		sum1 += d.i_ * (1 - d.j_) * tmp
		sum2 += d.i_ * (d.i_ + 1) * tmp
	}
	sum1 += n1_r3
	sum2 += n1_r3
	return math.sqrt(1000.0 * if97_r *  T * (sum2 - sq(sum1)/ttphitt))
}

//
pub fn r3_alphav(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in rt_r3 {
		tmp := d.n_ * d.i_ * ipow2(del, tau, d.i_, d.j_)
		sum1 += (1 - d.j_) * tmp
		sum2 += (d.i_ + 1) * tmp
	}
	sum1 += n1_r3
	sum2 += n1_r3
	return sum1/sum2/T
}

pub fn r3_kt(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut sum := 0.0
	for d in rt_r3 {
		sum += d.n_ * d.i_ * (d.i_ + 1) * ipow2(del, tau, d.i_, d.j_)
	}
	sum += n1_r3
    return 1.0/sum/rho/if97_r/T*1000.0
}

//
pub fn r3_alphap(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut tphit := 0.0
	mut dphid := 0.0
	for d in rt_r3 {
		tmp := d.n_ * ipow2(del, tau, d.i_, d.j_)
		tphit += d.j_ * tmp
		dphid += d.i_ * tmp
	}
	dphid += n1_r3
	return (1.0 - tphit/dphid) / T
}

pub fn r3_betap(rho f64, T f64) f64{
	del, tau := r3_deltau(rho, T)
	mut ddphidd := 0.0
	mut dphid := 0.0
	for d in rt_r3 {
		tmp := d.n_ * d.i_ * ipow2(del, tau, d.i_, d.j_)
		dphid += tmp
		ddphidd += (d.i_ - 1) * tmp
	}
	ddphidd += -n1_r3
	dphid += n1_r3
	return rho*(2.0 + ddphidd/dphid)
}

//****************************************************************
// REGION 4 G(p,T) EQUATIONS, 
// 273.15 K < T < 647.096 K
//****************************************************************
// psat(T), IF97 Eq 30, saturated line
pub fn r4_psat(T f64) f64{
	// 273.15 ≤ T ≤ 647.096 
	var_theta := T + pt_r4[8] / (T + pt_r4[9])

	aa := (var_theta + pt_r4[0]) * var_theta + pt_r4[1]
	bb := (pt_r4[2] * var_theta + pt_r4[3]) * var_theta + pt_r4[4]
	cc := (pt_r4[5] * var_theta + pt_r4[6]) * var_theta + pt_r4[7]
	
	tmp := 2.0 * cc / (- bb + math.sqrt(sq(bb) - 4.0 * aa * cc))
	return pstar_region245 * sq(sq(tmp))
}

// tsat(p), IF97 Eq 31, saturated line, APPROXIMATION Method
pub fn r4_tsat(p f64) f64{
	// 0.00061121 ≤ P ≤ 22.064
	pi := p / pstar_region245
	beta := math.pow(pi, 0.25)

	ee := (beta + pt_r4[2]) * beta + pt_r4[5]
	ff := (pt_r4[0] * beta + pt_r4[3]) * beta + pt_r4[6]
	gg := (pt_r4[1] * beta + pt_r4[4]) * beta + pt_r4[7]
	dd := 2.0 * gg / (-ff - math.sqrt(sq(ff) - 4.0 * ee * gg))

	return 0.5 * (dd - pt_r4[9] - math.sqrt(- 4.0 * pt_r4[8] + sq(dd + pt_r4[9]))) * tstar_region4
}

//------------------------------------------------------------------------------
// density
fn r4_rhof(T f64) f64{
	tau := 1.0 - T / if97_tcrit

	taui_3 := math.pow(tau, 1.0/3)

	tau_2_3 := sq(taui_3)
	tau_5_3 := tau * tau_2_3
	taui6_3 := sq(tau_5_3) * tau_5_3 * taui_3
	tau_43_3 := sq(taui6_3) * sq(tau_5_3) * taui_3
	taui10_3 := sq(tau_43_3) * taui6_3 * tau_5_3 * tau

	delta := 1 + t1_r4[0]*taui_3 + t1_r4[1]*tau_2_3 + t1_r4[2]*tau_5_3	+ t1_r4[3]*taui6_3	+ t1_r4[4]*tau_43_3	+ t1_r4[5]*taui10_3

	return delta * if97_rhocrit

}

fn r4_rhog(T f64) f64{
	tau := 1.0 - T / if97_tcrit

	taui_6 := math.pow(tau,1.0/6)

	tau_2_6 := sq(taui_6)
	tau_4_6 := sq(tau_2_6)
	tau_8_6 := sq(tau_4_6)
	taui6_6 := sq(tau_8_6)
	taui8_6 := taui6_6 * tau_2_6
	tau_37_6 := sq(taui8_6) * taui_6
	tau_71_6 := tau_37_6 * taui8_6 * taui6_6

	ln_delta := t2_r4[0]*tau_2_6 + t2_r4[1]*tau_4_6	+ t2_r4[2]*tau_8_6 + t2_r4[3]*taui8_6 + t2_r4[4]*tau_37_6 + t2_r4[5]*tau_71_6

	return math.exp(ln_delta) * if97_rhocrit
}

//****************************************************************
// REGION 5 G(p,T) EQUATIONS , IF97-Rev, Eq 32-34
// 1073.15 K < T < 2273.15 K, 0 < p < 50 MPa
//****************************************************************
// Internal function
fn r5_pitau(p f64, t f64) (f64, f64) {
	return p / pstar_region245, tstar_region5/ t
}

// Public function
pub fn r5_v(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum := 0.0
	for d in pt_r5{
		if d.i_ == 0 {continue}
		sum += d.n_ * d.i_ * ipow2(pi, tau, d.i_, d.j_)
	}
	return if97_r *  T * (sum + 1.0) / p / 1000.0
}

pub fn r5_u(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum := 0.0
	for d in pt_r5{
		sum += d.n_ * (d.j_ - d.i_) * ipow2(pi, tau, d.i_, d.j_) 
	}
	return if97_r * T * (sum-1.0)
}

pub fn r5_s(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum := 0.0
	for d in pt_r5{
		sum += d.n_ * (d.j_ - 1) * ipow2(pi, tau, d.i_, d.j_) 
	}
	return if97_r * (sum-math.log(pi))
}

pub fn r5_h(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum := 0.0
	for d in pt_r5{
		sum += d.n_ * d.j_ * ipow2(pi, tau, d.i_, d.j_) 
	}
	return if97_r * T * sum 
}

pub fn r5_cp(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum := 0.0
	for d in pt_r5{
		sum += d.n_ * d.j_ * (1 - d.j_) * ipow2(pi, tau, d.i_, d.j_) 
	}
	return if97_r * sum
}

pub fn r5_cv(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum1 := 0.0
	mut sum2 := 0.0
	mut sum3 := 0.0
	for d in pt_r5{
		tmp := d.n_ * ipow2(pi, tau, d.i_, d.j_)
		sum1 += d.j_ * (1 - d.j_) * tmp
		if d.i_ != 0 {
			sum2 += d.i_ * (1 - d.j_) * tmp 
			sum3 += d.i_ * (1 - d.i_) * tmp 
		}
	}
	return if97_r * (sum1-sq(1.0+sum2)/(1.0+sum3))
}

pub fn r5_w(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum1 := 0.0
	mut sum2 := 0.0
	mut sum3 := 0.0
	mut sum4 := 0.0
	for d in pt_r5{
		tmp := d.n_ * ipow2(pi, tau, d.i_, d.j_)
		sum1 += d.j_ * (d.j_ - 1) * tmp
		if d.i_ != 0 {
			tmp2 := d.i_ * tmp
			sum4 += tmp2
			sum2 += (1 - d.j_) * tmp2 
			sum3 += (1 - d.i_) * tmp2 
		}
	}
	return math.sqrt(1000.0*if97_r * T/(1.0+sum3+sq(1.0+sum2)/sum1))*(1.0+sum4)
}

pub fn r5_g(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum := 0.0
	for d in pt_r5{
		sum += d.n_ * ipow2(pi, tau, d.i_, d.j_)
	}
	return if97_r *  T * (sum + math.log(pi))
}

pub fn r5_a(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum := 0.0
	for d in pt_r5{
		sum += d.n_ * (1 - d.i_) * ipow2(pi, tau, d.i_, d.j_)
	}
	sum += (math.log(pi) - 1.0)
	return if97_r *  T * sum
}

pub fn r5_alphav(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in pt_r5{
		if d.i_ == 0 { continue }
		tmp := d.n_ * d.i_ *  ipow2(pi, tau, d.i_, d.j_)
		sum2 += tmp 
		sum1 += (1 - d.j_) * tmp
	}	
	return (1.0 + sum1) / (1.0 + sum2) / T
}

pub fn r5_kt(p f64, T f64) f64{
	pi, tau := r5_pitau(p, T)
	mut sum1 := 0.0
	mut sum2 := 0.0
	for d in pt_r5{
		if d.i_ == 0 { continue }
		tmp := d.n_ * d.i_ *  ipow2(pi, tau, d.i_, d.j_)
		sum2 += tmp 
		sum1 += (1 - d.i_) * tmp
	}	
	return (1.0 + sum1) / (1.0 + sum2) / p
}
