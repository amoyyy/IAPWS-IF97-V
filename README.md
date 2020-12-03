# IAPWS V Package

V module that provides IAPWS-IF97, IAPWS-IF95 support and some helper functions, WORK IN PROGRESS

References: http://iapws.org/release.html

Install:
        
	v install amoyyy.iapws

Usage:
	
	import amoyyy.iapws.if97
	t := 700.0	// K
	p := 30.0	// Mpa
	
	println(if97.r2_u(p, t))
