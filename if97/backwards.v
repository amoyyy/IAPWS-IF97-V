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
 * File: backwards.v
 * Author: YU YANG
 * Created Time: 2020-05-07
 * Version: 0.0.1
 *
 **********************************************************************/
 module if97
 import math

//****************************************************************
// REGION 1 backwards
// 273.15 K < T < 623.15 K, ps(T) < p < 100 MPa.
//****************************************************************
// T=f(P,h), IF97-Rev, Eq 11
pub fn bw1_t_ph(p f64, h f64) f64{
	pi := p/1.0
	yita_ := h/2500.0+1.0

	mut t := 0.0
	for d in tph_bw1{
		t += d.n_ * ipow2(pi, yita_, d.i_, d.j_)
	}

	return t
}

// T=f(P,s), IF97-Rev, Eq 13
pub fn bw1_t_ps(p f64, h f64) f64{
	pi := p/1.0
	sigma_ := h/1.0+2.0

	mut t := 0.0
	for d in tps_bw1{
		t += d.n_ * ipow2(pi, sigma_, d.i_, d.j_)
	}

	return t
}

// p=f(h,s), Supp-PHS12-2014.pdf, Eq 1
pub fn bw1_p_hs(h f64, s f64) f64{
	yita_ := h/3400.0+0.05
	sigma_ := s/7.6+0.05

	mut p := 0.0
	for d in phs_bw1{
		p += d.n_ * ipow2(yita_, sigma_, d.i_, d.j_)
	}

	return 100.0*p
}

// t=f(h,s), Supp-PHS12-2014.pdf, Section 7.1
pub fn bw1_t_hs(h f64, s f64) f64{
	p := bw1_p_hs(h, s)
	t := bw1_t_ph(p, h)
	return t
}

//****************************************************************
// REGION 2 backwards
// 273.15 K < T < 623.15 K, 0 < p < ps(T).
// 623.15 K < T < 863.15 K, 0 < p < b23_p(T).
// 863.15 K < T < 1073.15 K, 0 < p < 100 MPa.
//****************************************************************
// T=f(P,h), IF97-Rev, Eq 22
fn bw2a_t_ph(p f64, h f64) f64{
	pi := p/1.0
	sigma_ := h/2000.0-2.1
	mut t := 0.0
	for d in tph_bw2a{
		t += d.n_ * ipow2(pi, sigma_, d.i_, d.j_)
	}
	return t
}

// T=f(P,h), IF97-Rev, Eq 23
fn bw2b_t_ph(p f64, h f64) f64{
	pi_ := p/1.0 - 2.0
	yita_ := h/2000.0 - 2.6
	mut t := 0.0
	for d in tph_bw2b{
		t += d.n_ * ipow2(pi_, yita_, d.i_, d.j_)
	}
    return t
}

// T=f(P,h), IF97-Rev, Eq 24
fn bw2c_t_ph(p f64, h f64) f64{
	pi_ := p/1.0 + 25.0
	yita_ := h/2000.0 - 1.8
	mut t := 0.0
	for d in tph_bw2c{
		t += d.n_ * ipow2(pi_, yita_, d.i_, d.j_)
	}
    return t
}

// T=f(P,s), IF97-Rev, Eq 25
fn bw2a_t_ps(p f64, s f64) f64{
  	pi_ := p/1.0
	sigma_ := s/2.0 - 2.0
	mut t := 0.0
	for d in tps_bw2a{
		t += d.n_ * math.pow(pi_, d.i_)*ipow(sigma_, d.j_)
	}
    return t
}

// T=f(P,s), IF97-Rev, Eq 26
fn bw2b_t_ps(p f64, s f64) f64{
	pi_ := p/1.0
	sigma_ := 10.0 - s/0.7853
	mut t := 0.0
	for d in tps_bw2b{
		t += d.n_ * ipow2(pi_, sigma_, d.i_, d.j_)
	}
    return t
}

// T=f(P,s), IF97-Rev, Eq 27
fn bw2c_t_ps(p f64, s f64) f64{
	pi_ := p/1.0
	sigma_ := 2.0 - s/2.9251
	mut t := 0.0
	for d in tps_bw2c{
		t += d.n_ * ipow2(pi_, sigma_, d.i_, d.j_)
	}
    return t
}

// p=f(h,s), Supp-PHS12-2014.pdf, Eq 3
fn bw2a_p_hs(h f64, s f64) f64{
	yita_ := h/4200.0 - 0.5
	sigma_ := s/12.0 - 1.2
	mut p := 0.0
	for d in phs_bw2a{
		p += d.n_ * ipow2(yita_, sigma_, d.i_, d.j_)
	}
    return 4.0*sq(sq(p))
}

// p=f(h,s), Supp-PHS12-2014.pdf, Eq 4
fn bw2b_p_hs(h f64, s f64) f64{
	yita_ := h/4100.0 - 0.6
	sigma_ := s/7.9 - 1.01
	mut p := 0.0
	for d in phs_bw2b{
		p += d.n_ * ipow2(yita_, sigma_, d.i_, d.j_)
	}
    return 100.0*sq(sq(p))
}

// p=f(h,s), Supp-PHS12-2014.pdf, Eq 5
pub fn bw2c_p_hs(h f64, s f64) f64{
	yita_ := h/3500.0 - 0.7
	sigma_ := s/5.9 - 1.1
	mut p := 0.0
	for d in phs_bw2c{
		p += d.n_ * ipow2(yita_, sigma_, d.i_, d.j_)
	}
    return 100.0*sq(sq(p))
}

//-------------------------------------------------------------------
// T=f(P,h), Backward equation for region 2
pub fn bw2_t_ph(p f64, h f64) f64{
	mut t := 0.0
    if p <= 4.0{
        t = bw2a_t_ph(p, h)
	}else if p <= 6.546699678{ //T = 554.485 K and ps = 6.54670 MPa
        t = bw2b_t_ph(p, h)
	}else{
        if h >= b2bc_h_p(p){
        	t = bw2b_t_ph(p, h)
		}else{
        	t = bw2c_t_ph(p, h)
		}
	}

    if p <= 22.064 {
		tsat := r4_tsat(p)
		if t < tsat {return tsat}
	}
    return t
}

// T=f(P,s)
pub fn bw2_t_ps(p f64, s f64) f64{
	mut t := 0.0
    if p <= 4.0 {
		t = bw2a_t_ps(p,s)
	}else if s >= 5.85{
		t = bw2b_t_ps(p,s)
	}else{
		t = bw2c_t_ps(p,s)
	}

    if p <= 22.064 {
		tsat := r4_tsat(p)
		if t < tsat {return tsat}
	}
    return t
}

// p=f(h,s), for region2
pub fn bw2_p_hs(h f64, s f64) f64{
	s_2bc := 5.85 				// 5.85 kJ/kgK
    h_2ab := b2ab_h_s(s_2bc)		
    if h <= h_2ab {
		return bw2a_p_hs(h ,s)
	}else if s >= s_2bc{
		return bw2b_p_hs(h, s)
	}else{
		return bw2c_p_hs(h, s)
	}
}

// t=f(h,s), for region2
pub fn bw2_t_hs(h f64, s f64) f64{
	p := bw2_p_hs(h, s)		
	t := bw2_t_ph(p,h)
    return t
}

//****************************************************************
// REGION 3 backwards
//****************************************************************
// P=f(h,s), Supp-phs3-2014.pdf. Eq 1
fn bw3a_p_hs(h f64, s f64) f64{
	yita_ := h/2300.0 - 1.01
	sigma_ := s/4.4 - 0.750
	mut p := 0.0
	for d in phs_bw3a{
		p += d.n_ * ipow2(yita_, sigma_, d.i_, d.j_)
	}
    return 99.0*p
}

// P=f(h,s), Supp-phs3-2014.pdf. Eq 2
fn bw3b_p_hs(h f64, s f64) f64{
	yita_ := h/2800.0 - 0.681
	sigma_ := s/5.3 - 0.792
	mut p := 0.0
	for d in phs_bw3b{
		p += d.n_ * ipow2(yita_, sigma_, d.i_, d.j_)
	}
    return 16.6/p
}

// t=f(P,h), Supp-Tv(ph,ps)3-2014.pdf, Eq 2
fn bw3a_t_ph(p f64, h f64) f64{
    pi_ := p/100.0 + 0.240
	yita_ := h/2300.0 - 0.615
	mut t := 0.0
	for d in vph_bw3a{
		t += d.n_ * ipow2(pi_, yita_, d.i_, d.j_)
	}
	return 760.0*t
}

// t=f(P,h), Supp-Tv(ph,ps)3-2014.pdf, Eq 3
fn bw3b_t_ph(p f64, h f64) f64{
	pi_ := p/100.0 + 0.298
	yita_ := h/2800.0 - 0.720
	mut t := 0.0
	for d in vph_bw3a{
		t += d.n_ * ipow2(pi_, yita_, d.i_, d.j_)
	}
	return 860.0*t
}

// v=f(P,h), Supp-Tv(ph,ps)3-2014.pdf, Eq 4
fn bw3a_v_ph(p f64, h f64) f64{
    pi_ := p/100.0 + 0.128
	yita_ := h/2100.0 - 0.727
	mut v := 0.0
	for d in vph_bw3a{
		v += d.n_ * ipow2(pi_, yita_, d.i_, d.j_)
	}
	return 0.0028*v
}

// v=f(P,h), Supp-Tv(ph,ps)3-2014.pdf, Eq 5
fn bw3b_v_ph(p f64, h f64) f64{
	pi_ := p/100.0 + 0.0661
	yita_ := h/2800.0 - 0.720
	mut v := 0.0
	for d in vph_bw3b{
		v += d.n_ * ipow2(pi_, yita_, d.i_, d.j_)
	}
	return 0.0088*v
}

// T=f(P,s), Supp-Tv(ph,ps)3-2014.pdf, Eq 6
fn bw3a_t_ps(p f64, s f64) f64{
	pi_ := p/100.0 + 0.240
	sigma_ := s/4.4 - 0.703
	mut t := 0.0
	for d in tps_bw3a{
		t += d.n_ * ipow2(pi_, sigma_, d.i_, d.j_)
	}
	return 760.0*t
}

// T=f(P,s), Supp-Tv(ph,ps)3-2014.pdf, Eq 7
fn bw3b_t_ps(p f64, s f64) f64{
	pi_ := p/100.0 + 0.760
	sigma_ := s/5.3 - 0.818
	mut t := 0.0
	for d in tps_bw3b{
		t += d.n_ * ipow2(pi_, sigma_, d.i_, d.j_)
	}
	return 860.0*t
}

// v=f(P,s), Supp-Tv(ph,ps)3-2014.pdf, Eq 8
fn bw3a_v_ps(p f64, s f64) f64{
	pi_ := p/100.0 + 0.187
	sigma_ := s/4.4 - 0.755
	mut v := 0.0
	for d in vps_bw3a{
		v += d.n_ * ipow2(pi_, sigma_, d.i_, d.j_)
	}
	return 0.0028*v
}

// v=f(P,s), Supp-Tv(ph,ps)3-2014.pdf, Eq 9
fn bw3b_v_ps(p f64, s f64) f64{
	pi_ := p/100.0 + 0.298
	sigma_ := s/5.3 - 0.816
	mut v := 0.0
	for d in vps_bw3b{
		v += d.n_ * ipow2(pi_, sigma_, d.i_, d.j_)
	}
	return 0.0088*v
}

//-------------------------------------------------------------------
// P=f(h,s), region 3
fn bw3_p_hs(h f64, s f64) f64{
    sc := 4.41202148223476
    if s <= sc{
        return bw3a_p_hs(h, s)
	}else{
        return bw3b_p_hs(h, s)
	}
}

// v=f(P,h), region 3
fn bw3_v_ph(p f64, h f64) f64{
    h_3ab := b3ab_h_p(p)
    if h <= h_3ab {
        return bw3a_v_ph(p, h)
	}else{
        return bw3b_v_ph(p, h)
	}
}

// t=f(P,h), region 3
fn bw3_t_ph(p f64, h f64) f64{
	h_3ab := b3ab_h_p(p)
    if h <= h_3ab {
        return bw3a_t_ph(p, h)
	}else{
        return bw3b_t_ph(p, h)
	}
}

// t=f(P,s), region 3
fn bw3_t_ps(p f64, s f64) f64{
	sc := 4.41202148223476
    if s <= sc{
        return bw3a_t_ps(p, s)
	}else{
        return bw3b_t_ps(p, s)
	}
}

// v=f(P,s), region 3
fn bw3_v_ps(p f64, s f64) f64{
	sc := 4.41202148223476
    if s <= sc{
        return bw3a_v_ps(p, s)
	}else{
        return bw3b_v_ps(p, s)
	}
}

//****************************************************************
// SUB-REGION 3 backwards
//****************************************************************
// v=f(P,T,x), Supp-VPT3-2016.pdf, Eq. 4-5
pub fn bw3_v_pt(p f64, t f64, s string) f64{
	mut pi := p/r3s[s].p_ - r3s[s].a_
	mut tau := t/r3s[s].t_ - r3s[s].b_
	if r3s[s].c_r > 1 {pi = math.sqrt(pi)}					// for c = 0.5
	if r3s[s].d_r > 1 {tau = math.sqrt(math.sqrt(tau))}		// for d = 0.25

    mut v := 0.0
	for d in r3s[s].ijn_{
		v += d.n_ * ipow2(pi, tau, d.i_, d.j_)
	}

    if s == "n" {
		v = math.exp(v)
	}else if r3s[s].e_ > 1 {
		v = sq(sq(v))
	}
    return r3s[s].v_*v
}

//****************************************************************
// REGION 4 backwards
//****************************************************************
// T=f(h,s), Supp-phs3-2014.pdf. Eq 9
pub fn bw4_t_hs(h f64, s f64) f64{
	yita_ := h/2800.0 - 0.119
	sigma_ := s/9.2 - 1.07
	mut t := 0.0
	for d in ths_bw4{
		t += d.n_ * ipow2(yita_, sigma_, d.i_, d.j_)
	}
    return 550 * t
}
