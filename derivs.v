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
 * File: derivs.v
 * Author: YU YANG
 * Created Time: 2020-05-07
 * Version: 0.0.1
 *
 **********************************************************************/
 module if97
 import math

//****************************************************************
// REGION 1 derivatives
//****************************************************************

//****************************************************************
// REGION 2 derivatives
//****************************************************************
//----------------------------------------------------------------

//****************************************************************
// REGION 3 derivatives
//****************************************************************


//****************************************************************
// REGION 4 derivatives
//****************************************************************
fn r4_dpdt(t f64) f64{
	/* implicit differentiation */
	beta := math.pow(r4_psat(t)/pstar_region245, 0.25)
	theta := t + pt_r4[8] / (t + pt_r4[9])

	xbeta := (2.0*beta + pt_r4[2])*sq(theta) + (2.0*beta*pt_r4[0] + pt_r4[3])*theta + 2.0*pt_r4[1]*beta + pt_r4[4]
	xtheta := (2.0*theta + pt_r4[0])*sq(beta) + (2.0*pt_r4[2]*theta + pt_r4[3])*beta + 2.0*pt_r4[5]*theta + pt_r4[6] 

	dtheta_dt := 1.0 - pt_r4[8] / (t + pt_r4[9])
	dbeta_dtheta := -xtheta/xbeta
	dp_dbeta := 4.0*sq(beta)*beta*pstar_region245

	return dp_dbeta * dbeta_dtheta * dtheta_dt
}
