# IAPWS V Package

V module that provides IAPWS-IF97 and some helper functions, WORK IN PROGRESS

IAPWS-IF95 support will be provided later.

References: http://iapws.org/release.html

Install:
        
	v install amoyyy.iapws97

Usage:
	
	import amoyyy.iapws97.if97
	t := 700.0	// K
	p := 30.0	// Mpa
	
	println(if97.r2_u(p, t))
