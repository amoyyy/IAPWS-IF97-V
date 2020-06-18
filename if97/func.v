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
fn fabs(x f64) f64 {
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
