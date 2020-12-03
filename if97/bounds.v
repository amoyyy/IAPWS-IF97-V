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
 * File: bounds.v
 * Author: YU YANG
 * Created Time: 2020-05-07
 * Version: 0.0.1
 *
 **********************************************************************/
module if97
import math

//****************************************************************
// BOUNDARY: region 1/3 (saturated line)
//****************************************************************
// h=f(s), Supp-phs3-2014.pdf. Eq 7, Boundary between Regions 1 and 3
pub fn b13_h_s(s f64) f64{
	// s(100MPa,623.15K) â‰? s â‰? s'(623.15K)
    error_input := s < 3.397782955 || s > 3.77828134
	
	if error_input {
		//error("B13_H_S input out of bound")
		return 0.0
	}else{
		sigma_1 := s/3.8 - 0.884
		sigma_2 := s/3.8 - 0.864
		mut sum := 0.0
		for d in s_b13 {
			sum += d.n_ * ipow2(sigma_1, sigma_2, d.i_, d.j_)
		}
		return 1700.0 * sum
	}
}

//****************************************************************
// BOUNDARY: region 1/4 (saturated line)
//****************************************************************
// h=f(s), Supp-phs3-2014.pdf. Eq 3
pub fn b14_h_s(s f64) f64{
	// s'(273.15K) â‰? s â‰? s'(623.15K)
	error_input := s < -1.545495919e-4 || s > 3.77828134
	
	if error_input{
		//error("B14_H1_S input out of bound")
		return 0.0
	}else{
		sigma_1 := s/3.8 - 1.09
		sigma_2 := s/3.8 + 0.366e-4
		mut sum := 0.0
		for d in s_b14 {
			sum += d.n_ * ipow2(sigma_1, sigma_2, d.i_, d.j_)
		}
		return 1700.0 * sum
	}
}

//****************************************************************
// BOUNDARY: region 2/3
//****************************************************************
// p=f(T), IF97-Rev, Eq 5, boundary between Region 2 and 3
pub fn b23_p_t(T f64) f64{
	return ((tp_b23[2] * T + tp_b23[1]) * T + tp_b23[0]) * pstar_region245
}

// T=f(P), IF97-Rev, Eq 6, boundary between Region 2 and 3
pub fn b23_t_p(p f64) f64{
	pi := p / pstar_region245
	return tp_b23[3] + math.sqrt((pi - tp_b23[4])/tp_b23[2])
}

// T=f(h,s), Supp-phs3-2014.pdf. Eq 8, boundary between Region 2 and 3
pub fn b23_t_hs(h f64, s f64) f64{
	// 5.048096828 â‰? s â‰? 5.260578707
    // 2.563592004e3 â‰? h â‰? 2.812942061e3
    error_input := h < 2.563592004e3 || h > 2.812942061e3 || s < 5.048096828 || s > 5.260578707

	if error_input{
		//error("B23_T_HS input out of bound")
		return 0.0
	}else{
		yita_ := h/3000.0 - 0.727
		sigma_ := s/5.3 - 0.864
		mut sum := 0.0
		for d in hs_b23 {
			sum += d.n_ * ipow2(yita_, sigma_, d.i_, d.j_)
		}
		return 900.0*sum
	}
}

//****************************************************************
// BOUNDARY: region 2b/2c, 2a/2b
//****************************************************************
// P=f(h), IF97-Rev, Eq 20, boundary between Region 2b and 2c
pub fn b2bc_p_h(h f64) f64{
    return hp_b2bc[0] + (hp_b2bc[1] + hp_b2bc[2] * h) * h
}

// h=f(P), IF97-Rev, Eq 21, boundary between Region 2b and 2c
pub fn b2bc_h_p(p f64) f64{
    return hp_b2bc[3] + math.sqrt((p-hp_b2bc[4])/hp_b2bc[2])
}

// h=f(s), Supp-PHS12-2014.pdf, Eq 2, boundary between Region 2a and 2b
pub fn b2ab_h_s(s f64) f64{
    // smin = r2_s(4MPa, tsat(4Mpa))
    // smax = r2_s(4MPa, 1073.15K)
    error_input := s < 6.069709159519128 || s > 7.85234039987851

    if error_input{
        //error("B2AB_H_S input not valid")
		return 0.0
    }
    return ((hs_b2ab[3]*s+hs_b2ab[2])*s+hs_b2ab[1])*s+hs_b2ab[0]
}

//****************************************************************
// BOUNDARY: REGION 3, 3a/3b, 3o/3p, 3w/3x, 3e/3f, 3x/3y
//****************************************************************
// P=f(h), Supp-Tv(ph,ps)3-2014.pdf, Eq 10, (Region3 Saturated Line)
pub fn b3_psat_h(h f64) f64{
	// h'(623.15K) â‰? h â‰? h''(623.15K)
    error_input := h < 1.670858218e3 || h > 2.563592004e3

	if error_input{
		//error("B3_PSAT_H input out of bound")
		return 0.0
	}else{
		yita_1 := h/2600.0 - 1.02
		yita_2 := h/2600.0 - 0.608
		mut sum := 0.0
		for d in h_b3 {
			sum += d.n_ * ipow2(yita_1, yita_2, d.i_, d.j_)
		}
		return 22.0*sum
	}
}

// P=f(s), Supp-Tv(ph,ps)3-2014.pdf, Eq 11, (Saturated Line)
pub fn b3_psat_s(s f64) f64{
    // r1_s(623.15K) â‰? s â‰? r2_s(623.15K)
	error_input := s < 3.778281340 || s > 5.210887825

	if error_input{
		//error("B3_PSAT_S input out of bound")
		return 0.0
	}else{
		sigma_1 := s/5.2 - 1.03
		sigma_2 := s/5.2 - 0.699
		mut sum := 0.0
		for d in s_b3 {
			sum += d.n_ * ipow2(sigma_1, sigma_2, d.i_, d.j_)
		}
		return 22.0*sum
	}
}

// h=f(P), Supp-Tv(ph,ps)3-2014.pdf, Eq 1, boundary between Region 3a/3b
pub fn b3ab_h_p(p f64) f64{
	// if97_pcrit â‰? p â‰? 100.0Mpa
	return hp_b3ab[0]+(hp_b3ab[1]+(hp_b3ab[2]+hp_b3ab[3]*p)*p)*p
}
