module if95

const (
    t_crit = 647.096            // K
    t_triple = 273.16           // K
    h_triple = 0.611782         // J/kg
    rho_crit = 322.0            // kg/m3
    rconst = 0.46151805         // kJ/(kg*K)
)

/* Data struct used throughout IAPWS-IF95 */
struct DTN{
    d int
    t f64
    n f64
	dd f64
	tt f64
	dt f64
}

struct CDTN{
    c int
    d int
    t int
	n f64
	tt f64
	c2 f64
}

fn dtn(d int, t f64, n f64) DTN { 
	return DTN{d: d, t:t, n:n, dd:f64(d*(d-1)), tt:t*(t-1), dt:f64(d)*t} 
}

fn cdtn(c int, d int, t int, n f64) CDTN { 
	return CDTN{c:c, d: d, t:t, n:n, tt:f64(t*(t-1)), c2:f64(c*c)} 
}

const (
    // phi0
    phi0_n = [0., -8.3204464837497, 6.6832105275932, 3.00632, 0.012436, 0.97315, 1.27950, 0.96956, 0.24873]
    phi0_gamma = [0., 0., 0., 0., -1.28728967, -3.53734222, -7.74073708, -9.24437796, -27.5075105]
    // phi1
	/*
[0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6]
[1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6]
[-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1, 4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23, 23, 10, 50, 44, 46, 50]
[0.012533547935523, 7.895763472282800, -8.780320330356100, 0.318025093454180, -0.261455338593580, -0.007819975168798, 0.008808949310213, -0.668565723079650, 0.204338109509650, -0.000066212605040, -0.192327211560020, -0.257090430034380, 0.160748684862510, -0.040092828925807, 0.000000393434226, -0.000007594137709, 0.000562509793519, -0.000015608652257, 0.000000001153800, 0.000000365821651, -0.000000000001325, -0.000000000626396, -0.107936009089320, 0.017611491008752, 0.221322951675460, -0.402476697635280, 0.580833999857590, 0.004996914699081, -0.031358700712549, -0.743159297103410, 0.478073299154800, 0.020527940895948, -0.136364351103430, 0.014180634400617, 0.008332650488071, -0.029052336009585, 0.038615085574206, -0.020393486513704, -0.001655405006373, 0.001995557197954, 0.000158703083242, -0.000016388568343, 0.043613615723811, 0.034994005463765, -0.076788197844621, 0.022446277332006, -0.000062689710415, -0.000000000557111, -0.199057183544080, 0.317774973307380, -0.118411824259810]

	*/
    phir_1 = [  dtn(1, -0.5, 0.12533547935523E-1),
                dtn(1, 0.875, 0.78957634722828E1),
                dtn(1, 1.0, -0.87803203303561E1),
                dtn(2, 0.5, 0.31802509345418),
                dtn(2, 0.75, -0.26145533859358),
                dtn(3, 0.375, -0.78199751687981E-2),
                dtn(4, 1.0, 0.88089493102134E-2) 
			 ]
	phir_2 = [  cdtn(1, 1, 4, -0.66856572307965),
				cdtn(1, 1, 6, 0.20433810950965),
				cdtn(1, 1, 12, -0.66212605039687E-4),
				cdtn(1, 2, 1, -0.19232721156002),
				cdtn(1, 2, 5, -0.25709043003438),
				cdtn(1, 3, 4, 0.16074868486251),
				cdtn(1, 4, 2, -0.40092828925807E-1),
				cdtn(1, 4, 13, 0.39343422603254E-6),
				cdtn(1, 5, 9, -0.75941377088144E-5),
				cdtn(1, 7, 3, 0.56250979351888E-3),
				cdtn(1, 9, 4, -0.15608652257135E-4),
				cdtn(1, 10, 11, 0.11537996422951E-8),
				cdtn(1, 11, 4, 0.36582165144204E-6),
				cdtn(1, 13, 13, -0.13251180074668E-11),
				cdtn(1, 15, 1, -0.62639586912454E-9),
				cdtn(2, 1, 7, -0.10793600908932),
				cdtn(2, 2, 1, 0.17611491008752E-1),
				cdtn(2, 2, 9, 0.22132295167546),
				cdtn(2, 2, 10, -0.40247669763528),
				cdtn(2, 3, 10, 0.58083399985759),
				cdtn(2, 4, 3, 0.49969146990806E-2),
				cdtn(2, 4, 7, -0.31358700712549E-1),
				cdtn(2, 4, 10, -0.74315929710341),
				cdtn(2, 5, 10, 0.47807329915480),
				cdtn(2, 6, 6, 0.20527940895948E-1),
				cdtn(2, 6, 10, -0.13636435110343),
				cdtn(2, 7, 10, 0.14180634400617E-1),
				cdtn(2, 9, 1, 0.83326504880713E-2),
				cdtn(2, 9, 2, -0.29052336009585E-1),
				cdtn(2, 9, 3, 0.38615085574206E-1),
				cdtn(2, 9, 4, -0.20393486513704E-1),
				cdtn(2, 9, 8, -0.16554050063734E-2),
				cdtn(2, 10, 6, 0.19955571979541E-2),
				cdtn(2, 10, 9, 0.15870308324157E-3),
				cdtn(2, 12, 8, -0.16388568342530E-4),
				cdtn(3, 3, 16, 0.43613615723811E-1),
				cdtn(3, 4, 22, 0.34994005463765E-1),
				cdtn(3, 4, 23, -0.76788197844621E-1),
				cdtn(3, 5, 23, 0.22446277332006E-1),
				cdtn(4, 14, 10, -0.62689710414685E-4),
				cdtn(6, 3, 50, -0.55711118565645E-9),
				cdtn(6, 6, 44, -0.19905718354408),
				cdtn(6, 6, 46, 0.31777497330738),
				cdtn(6, 6, 50, -0.11841182425981),
			 ]
/*
	phir_3 = [
				(3, 0, -0.31306260323435E2, 20., 150., 1.21, 1.0),
				(3, 1, 0.31546140237781E2, 20., 150., 1.21, 1.0),
				(3, 4, -0.25213154341695E4, 20., 250., 1.25, 1.0)
	]
	phir_4 = [
				(3.5, 0.85, 0.2, -0.14874640856724, 28, 700, 0.32, 0.3),
			    (3.5, 0.95, 0.2, 0.31806110878444, 32, 800, 0.32, 0.3)
	]
*/
)

