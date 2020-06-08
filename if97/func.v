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
 * File: func.v
 * Author: YU YANG
 * Created Time: 2020-05-26
 * Version: 0.0.1
 *
 **********************************************************************/
 module if97

/* Number Functions*/
/* square */
fn sq(x f64) f64 {
	return x*x
}

/* cubic */
fn cube(x f64) f64 {
	return x*x*x
}

/* fabs */ 
pub fn fabs(x f64) f64 {
	//*(((int *) &x) + 1) &= 0x7fffffff;
	if x < 0 {return -x}
	return x
}

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

/* ipow:  public domain by Mark Stephen with suggestions by Keiichi Nakasato */
fn ipow2(x1_in f64, x2_in f64, n1_in int, n2_in int) f64{
	mut x1 := x1_in
	mut n1 := n1_in 
	mut x2 := x2_in
	mut n2 := n2_in
	mut t := 1.0

	if n1 < 0{
		n1 = -n1
		x1 = 1.0/x1
	}

	x1 *= x2
	n2 -= n1

	//if n2 == 0 {return 1.0}

	if n2 < 0{
		n2 = -n2
		x2 = 1.0/x2
	}

	for n1 > 0 {
		if n1&1==1 {t*=x1}
		n1 /= 2     /* KN prefers if (n/=2) x*=x; This avoids an    */
		x1 *= x1     /* unnecessary but benign multiplication on     */
	}     		   /* the last pass, but the comparison is always
					true _except_ on the last pass. */

	for n2 > 0 {
		if n2&1==1 {t*=x2}
		n2 /= 2   
		x2 *= x2   
	}  

	return t
}

/*
Rejected since issue #5082: cannot define methods on types from other modules
/* Polynomial ops  YU YANG*/
fn (coeff []f64) po_mul(vars f64) f64 {
	mut sum := 0.0
	for i:=0; i < coeff.len - 1; i++{
		sum += coeff[i]
		sum *= vars
	}
	sum += coeff[coeff.len-1]
	return sum
}

fn (coeff []f64) po_mul_reverse(vars f64) f64 {
	mut sum := 0.0
	for i:=coeff.len - 1; i >= 0; i--{
		sum *= vars
		sum += coeff[i]
	}
	return sum
}

/* Array ops  YU YANG*/
fn (mut a []f64) add(input []f64){
	//assert a.len == input.len
	for i:=0; i < a.len; i++{
		a[i] += input[i]
	}
}

fn (mut a []f64) sub(input []f64){
	assert a.len == input.len
	for i:=0; i < a.len; i++{
		a[i] -= input[i]
	}
}

fn (mut a []f64) mul(input f64){
	for i:=0; i < a.len; i++{
		a[i] *= input
	}
}




fn r1_cal_pow(p f64, T f64) ([]f64, []f64) {
	pi, tau := r1_pitau(p, T)
	pii := 7.1 - pi 
	mut ipow_ := []f64{len: 33}
	mut jpow_ := []f64{len: 42}
	ipow_[0] = 1.0
	ipow_[1] = pii
	ipow_[2] = sq(pii)
	ipow_[3] = ipow_[1]*ipow_[2]
	ipow_[4] = ipow_[2]*ipow_[2]
	ipow_[5] = ipow_[2]*ipow_[3]
	ipow_[8] = ipow_[4]*ipow_[4]
	ipow_[21] = ipow_[8]*ipow_[8]*ipow_[5]
	ipow_[23] = ipow_[2]*ipow_[21]
	ipow_[29] = ipow_[23]*ipow_[5]*pii
	ipow_[30] = ipow_[29]*pii
	ipow_[31] = ipow_[30]*pii
	ipow_[32] = ipow_[31]*pii

	mut jmul := 1.0
	for j:=0;j<42;j++{
		jpow_[j] = jmul 
		jmul *= tau - 1.222
	}
	return ipow_, jpow_
}
*/