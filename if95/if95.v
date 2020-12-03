/**********************************************************************
 *
 * Vector Space System Project / Material Library Module / IAPWS-IF95
 *
 * Copyright (C) 2020 CIAE.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 * File: if95.v
 * Author: YU YANG
 * Created Time: 2020-10-31
 * Version: 0.0.1
 *
 **********************************************************************/
module if95
import math

pub struct IF95{
mut:
	name 	string = "IAPWS-95"
	del_pow []f64 = []f64{init:1.0, len:15}
	tau_pow []f64 = []f64{init:1.0, len:50}
	del		f64
	tau		f64
	// results
	phi 	f64
	phi_d 	f64
	phi_dd 	f64
	phi_t 	f64
	phi_tt 	f64
	phi_dt 	f64
	// properties
pub mut:
	rho			f64
	temperature f64
	p 		f64
    u		f64
    s		f64
    h		f64 
    cv 		f64
	cp		f64
	w		f64
}

pub fn (mut if95 IF95) set_state(rho f64, t f64) {
	if95.rho = rho
	if95.temperature = t
	del, tau := deltau(rho, t)
	if95.del = del
	if95.tau = tau
	if95.cal_init()
	if95.cal()
}

fn (mut if95 IF95) cal_init() {
	// init pow table
	mut mul_del := 1.0
	mut mul_tau := 1.0
	for idx := 1; idx <= 14; idx ++ {
		mul_del *= if95.del
		if95.del_pow[idx] = mul_del
	}
	for idx := 1; idx <= 12; idx ++ {
		mul_tau *= if95.tau
		if95.tau_pow[idx] = mul_tau
	}
	if95.tau_pow[15] = if95.tau_pow[7]*if95.tau_pow[8]
	if95.tau_pow[21] = if95.tau_pow[10] * if95.tau_pow[11]
	if95.tau_pow[22] = sq(if95.tau_pow[11])
	if95.tau_pow[43] = if95.tau_pow[21]*if95.tau_pow[22]
	if95.tau_pow[45] = if95.tau_pow[2]*if95.tau_pow[43]
	if95.tau_pow[49] = if95.tau_pow[4] * if95.tau_pow[45]
	// init result data
	if95.phi = 0.0
	if95.phi_d = 0.0
	if95.phi_dd = 0.0
	if95.phi_t = 0.0 
	if95.phi_tt = 0.0 
	if95.phi_dt = 0.0
}

[inline]
fn (mut if95 IF95) dt() (f64, f64) {
	return if95.del, if95.tau
}

[inline]
fn deltau(rho f64, t f64) (f64, f64) {
	return rho / rho_crit, t_crit / t
}

//*********************************************************************
// number functions
//*********************************************************************
/* square */
[inline]
fn sq(x f64) f64 {
	return x*x
}

[inline]
/* cubic */
fn cube(x f64) f64 {
	return x*x*x
}

[inline]
/* ipow:  public domain by Mark Stephen with suggestions by Keiichi Nakasato */
fn ipow(x_in f64, n_in int) f64{
	mut x := x_in
	mut n := n_in
	mut t := 1.0

	if n==0 {return 1.0}    /* At the top. x^0 = 1 */

	if n < 0{
		n = -n
		x = 1.0/x    /* error if x == 0. Good                        */
	}                /* ZTC/SC returns inf, which is even better     */

	if x == 0.0 {return 0.0}

	for n > 0 {
		if n&1==1 {t*=x}
		n /= 2     /* KN prefers if (n/=2) x*=x; This avoids an    */
		x *= x     /* unnecessary but benign multiplication on     */
	}     		   /* the last pass, but the comparison is always
					true _except_ on the last pass. */

	return t
}

//*********************************************************************
// phi_0 and derivatives
//*********************************************************************
[direct_array_access]
fn (mut if95 IF95) phi0() {
	del, tau := if95.dt()
    mut phi0 := -8.3204464837497 + math.log(del) + 6.6832105275932*tau + 3.00632*math.log(tau)
    mut phi0_t := 6.6832105275932 + 3.00632 / tau
    mut phi0_tt := - 3.00632 / sq(tau)
    for ii:=4; ii<=8; ii++ {
        tmp := math.exp(phi0_gamma[ii]*tau)
        tmp2 := 1.0 - tmp
		tmp3 := - phi0_n[ii] * phi0_gamma[ii] * tmp / tmp2
        phi0 += phi0_n[ii] * math.log(tmp2)
        phi0_t += tmp3
        phi0_tt += tmp3 * phi0_gamma[ii] / tmp2
    }
	if95.phi += phi0
    if95.phi_d += 1.0/del
    if95.phi_dd += -1.0/sq(del)
    if95.phi_t += phi0_t
    if95.phi_tt += phi0_tt
	//if95.phi_dt += 0.0
}

//*********************************************************************
// phi_r and derivatives
//*********************************************************************
[direct_array_access]
fn (mut if95 IF95) phir() {
	mut phir := 0.0
	mut phir_d := 0.0
	mut phir_dd := 0.0
	mut phir_t := 0.0
	mut phir_tt := 0.0
	mut phir_dt := 0.0

	del, tau := if95.dt()

	// part 1 
    for dtn in phir_1 {
		d := dtn.d
		t := dtn.t
		deld := if95.del_pow[d-1]

        tmp := dtn.n * deld * math.pow(tau, t-1.0)
        
		phir += tmp
        phir_d += tmp * f64(d)
        phir_dd += tmp * dtn.dd
        phir_t += tmp * t
        phir_tt += tmp * dtn.tt
        phir_dt += tmp * dtn.dt
    }
	// part 2
	for cdtn in phir_2 {
		c := cdtn.c
		d := cdtn.d
		t := cdtn.t
        delc := if95.del_pow[c]
		deld := if95.del_pow[d-1]
		taut := if95.tau_pow[t-1]
        para_d := f64(d) - f64(c) * delc
		para_dd := para_d * ( para_d - 1.0 ) - cdtn.c2 * delc

        tmp := cdtn.n * deld * taut * math.exp(-delc)

        phir += tmp
        phir_d += tmp * para_d
        phir_dd += tmp * para_dd
        phir_t += tmp * f64(t)
        phir_tt += tmp * cdtn.tt
        phir_dt += tmp * para_d * f64(t)
    }

	d_t := del/tau
	t_d := 1.0/d_t
	dt  := tau*del

	if95.phi += phir * dt
	if95.phi_d += phir_d * tau
	if95.phi_dd += phir_dd * t_d
	if95.phi_t += phir_t * del
	if95.phi_tt += phir_tt * d_t
	if95.phi_dt += phir_dt 

	// part3
	del_s1 := del - 1.0
	del_s2 := sq(del_s1)
    tau_s1 := tau - 1.0
    tau_s2 := sq(tau_s1)
    del_d1 := 1.0/del
    tau_d1 := 1.0/tau
    tau_d2 := sq(tau_d1)

    dpow3_d20t150 := if95.del_pow[3] * math.exp(-20.0*del_s2 - 150.0*sq(tau-1.21))
    tmp1 := -0.31306260323435E2 * dpow3_d20t150
    tmp2 := 0.31546140237781E2 * tau * dpow3_d20t150
    tmp3 := -0.25213154341695E4 * if95.del_pow[3] * if95.tau_pow[4] * math.exp(-20.0*del_s2 - 250.0*sq(tau-1.25)) 

    para_t1 := -300.0*tau + 363.0
    para_t2 := para_t1 + tau_d1
    para_t3 := -500.0*(tau-1.25) + 4.0*tau_d1

    phir = tmp1 + tmp2 + tmp3
    phir_t = tmp1 * para_t1 + tmp2 * para_t2 + tmp3 * para_t3
    phir_tt = tmp1 * (sq(para_t1) - 300.0) + tmp2 * (sq(para_t2) - 300.0 - tau_d2) + tmp3 * (sq(para_t3) - 500.0  - 4.0*tau_d2)
	
    para_d := 3.0*del_d1-40.0*del+40.0
    para_dd := -280.0+1600*del_s2+240.0*del_d1+6.0*sq(del_d1)

	if95.phi += phir
    if95.phi_d += phir * para_d
    if95.phi_dd += phir * para_dd
    if95.phi_t += phir_t 
    if95.phi_tt += phir_tt
    if95.phi_dt += phir_t * para_d

	// part4
	del_s2_7_6 := math.pow(del_s2, 7.0/6.0)
	del_s2_2_3 := math.pow(del_s2, 2.0/3.0)

    theta := 0.32 * del_s2 * del_s2_2_3 - tau_s1
    delta := sq(theta) + 0.2*cube(del_s2_7_6)
    delta_d := 32.0/15.0*theta*del_s2_7_6+1.4*cube(del_s2)
    delta_dd := del_s2_2_3*(del_s2_2_3*(8.4*del_s2_7_6 + 512.0/225.0*del_s2)+224.0/45.0*theta)

    tmp4 := -0.14874640856724 * math.pow(delta, 0.85) * del * math.exp(-28.0*del_s2-700.0*tau_s2)
    tmp5 := 0.31806110878444 * math.pow(delta, 0.95) * del * math.exp(-32.0*del_s2-800.0*tau_s2)

	c_phir := 56.0*tmp4 + 64.0*tmp5
	c2_phir := 3136.0*tmp4 + 4096.0*tmp5
	b_phir := (0.85*tmp4 + 0.95*tmp5)/delta
	bc_phir := (95.2*tmp4 + 121.6*tmp5)/delta
	bb_phir := (0.1275*tmp4 + 0.0475*tmp5)/sq(delta)
	bd_phir := (1190*tmp4 + 1520.0*tmp5)/delta
	d_phir := 1400.0*tmp4 + 1600.0*tmp5
	dd_phir := 1960000*tmp4 + 2560000*tmp5
	cd_phir := 78400.0*tmp4 + 102400.0*tmp5

	tmp6 := 2*b_phir*del_d1-bc_phir*del_s1
    
	phir = tmp4 + tmp5
    phir_d = phir*del_d1 - c_phir*del_s1 + b_phir*delta_d
    phir_t = -2.0*b_phir*theta - d_phir*tau_s1
    phir_dd = c_phir*(2.0*del_d1-3.0) + c2_phir*del_s2 + tmp6*delta_d + b_phir*delta_dd - bb_phir*sq(delta_d)
    phir_tt = 4.0*bd_phir*theta*tau_s1 + 2*b_phir - 4*bb_phir*sq(theta) + dd_phir*sq(tau_s1) - d_phir
	phir_dt = (-d_phir*del_d1 + cd_phir*del_s1 - bd_phir*delta_d)*tau_s1 - tmp6*theta - 32.0/15.0*b_phir*del_s2_7_6 + 2.0*bb_phir*theta*delta_d

	if95.phi += phir
    if95.phi_d += phir_d 
    if95.phi_dd += phir_dd
    if95.phi_t += phir_t 
    if95.phi_tt += phir_tt
    if95.phi_dt += phir_dt
}

//*************************************************************************************
// Property Calculation
//*************************************************************************************
fn (mut if95 IF95) cal(){
	if95.phi0()
	if95.phir()

	del, tau := if95.dt()
	d_t := del/tau

	u := if95.phi
	ut := if95.phi_t
	dt_udt := if95.phi_dt * del * tau
	t2_utt := if95.phi_tt * sq(tau)
	d_ud := if95.phi_d * del
	d2_udd := if95.phi_dd * sq(del)

    if95.p = rconst*rho_crit*t_crit*d_ud*d_t
	
    if95.u = rconst*t_crit*ut
    if95.s = rconst*(tau*ut-u)
    if95.h = rconst*t_crit*(ut+d_ud/tau)
    if95.cv = -rconst*t2_utt

    tmp1 := sq(d_ud-dt_udt)
    tmp2 := 2*d_ud+d2_udd

    if95.cp = if95.cv + rconst*tmp1/tmp2
    if95.w = math.sqrt(1000.0*rconst*t_crit/tau*(tmp2 - tmp1/t2_utt))
}

//*************************************************************************************
// Property Display
//*************************************************************************************
pub fn (if95 IF95) print_property(){
    print("    $if95.temperature    $if95.rho    ")
    println("$if95.p    $if95.cv    $if95.w    $if95.s    $if95.h")
}

pub fn (if95 &IF95) rho() f64{
	return if95.rho
}

pub fn (if95 &IF95) temperature() f64{
	return if95.temperature
}

pub fn (if95 &IF95) h() f64{
	return if95.h
}

pub fn (if95 &IF95) w() f64{
	return if95.w
}

pub fn (if95 &IF95) cp() f64{
	return if95.h
}

pub fn (if95 &IF95) s() f64{
	return if95.s
}

pub fn (if95 &IF95) p() f64{
	return if95.p
}
/*
pub fn test_if95(){
    // validation
	mut if95 := IF95{}
    print("   T [K]  ")
    print(" rho [kg/m3] ")
    print("    p [kpa]   ")
    print("cv [kJ/(kg*K)]")
    print("    w [m/s]   ")
    print(" s [kJ/(kg*K)]")
    println("   h [kJ/kg]  ")   

    test_t := [300., 500., 647., 900.]
    test_rho := [[996.556, 1005.308, 1188.202], 
	        [0.435, 4.532, 838.025, 1084.564], 
            [358.0],
            [0.241, 52.615, 870.769]]

    for mm:=0; mm<1; mm++{
		for ii:=0; ii<test_t.len; ii++{
			for jj:=0; jj<test_rho[ii].len; jj++{
				if95.set_state(test_rho[ii][jj], test_t[ii])
				if95.print_property()
			}
		}
    }
}
*/