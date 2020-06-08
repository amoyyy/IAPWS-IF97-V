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

const (
	tph_bw1 = [	ijn(0, 0, -2.3872489924521E+02),
				ijn(0, 1, 4.0421188637945E+02),
				ijn(0, 2, 1.1349746881718E+02),
				ijn(0, 6, -5.8457616048039E+00),
				ijn(0, 22, -1.5285482413140E-04),
				ijn(0, 32, -1.0866707695377E-06),
				ijn(1, 0, -1.3391744872602E+01),
				ijn(1, 1, 4.3211039183559E+01),
				ijn(1, 2, -5.4010067170506E+01),
				ijn(1, 3, 3.0535892203916E+01),
				ijn(1, 4, -6.5964749423638E+00),
				ijn(1, 10, 9.3965400878363E-03),
				ijn(1, 32, 1.1573647505340E-07),
				ijn(2, 10, -2.5858641282073E-05),
				ijn(2, 32, -4.0644363084799E-09),
				ijn(3, 10, 6.6456186191635E-08),
				ijn(3, 32, 8.0670734103027E-11),
				ijn(4, 32, -9.3477771213947E-13),
				ijn(5, 32, 5.8265442020601E-15),
				ijn(6, 32, -1.5020185953503E-17)]
	tps_bw1 = [	ijn(0, 0, 1.7478268058307E+02),
				ijn(0, 1, 3.4806930892873E+01),
				ijn(0, 2, 6.5292584978455E+00),
				ijn(0, 3, 3.3039981775489E-01),
				ijn(0, 11, -1.9281382923196E-07),
				ijn(0, 31, -2.4909197244573E-23),
				ijn(1, 0, -2.6107636489332E-01),
				ijn(1, 1, 2.2592965981586E-01),
				ijn(1, 2, -6.4256463395226E-02),
				ijn(1, 3, 7.8876289270526E-03),
				ijn(1, 12, 3.5672110607366E-10),
				ijn(1, 31, 1.7332496994895E-24),
				ijn(2, 0, 5.6608900654837E-04),
				ijn(2, 1, -3.2635483139717E-04),
				ijn(2, 2, 4.4778286690632E-05),
				ijn(2, 9, -5.1322156908507E-10),
				ijn(2, 31, -4.2522657042207E-26),
				ijn(3, 10, 2.6400441360689E-13),
				ijn(3, 32, 7.8124600459723E-29),
				ijn(4, 32, -3.0732199903668E-31)]
	phs_bw1 = [	ijn(0, 0, -6.91997014660582E-01),
				ijn(0, 1, -1.83612548787560E+01),
				ijn(0, 2, -9.28332409297335E+00),
				ijn(0, 4, 6.59639569909906E+01),
				ijn(0, 5, -1.62060388912024E+01),
				ijn(0, 6, 4.50620017338667E+02),
				ijn(0, 8, 8.54680678224170E+02),
				ijn(0, 14, 6.07523214001162E+03),
				ijn(1, 0, 3.26487682621856E+01),
				ijn(1, 1, -2.69408844582931E+01),
				ijn(1, 4, -3.19947848334300E+02),
				ijn(1, 6, -9.28354307043320E+02),
				ijn(2, 0, 3.03634537455249E+01),
				ijn(2, 1, -6.50540422444146E+01),
				ijn(2, 10, -4.30991316516130E+03),
				ijn(3, 4, -7.47512324096068E+02),
				ijn(4, 1, 7.30000345529245E+02),
				ijn(4, 4, 1.14284032569021E+03),
				ijn(5, 0, -4.36407041874559E+02)]
	tph_bw2a = [ijn(0, 0, 1089.8952318288),
				ijn(0, 1, 849.51654495535),
				ijn(0, 2, -107.81748091826),
				ijn(0, 3, 33.153654801263),
				ijn(0, 7, -7.4232016790248),
				ijn(0, 20, 11.765048724356),
				ijn(1, 0, 1.844574935579),
				ijn(1, 1, -4.1792700549624),
				ijn(1, 2, 6.2478196935812),
				ijn(1, 3, -17.344563108114),
				ijn(1, 7, -200.58176862096),
				ijn(1, 9, 271.96065473796),
				ijn(1, 11, -455.11318285818),
				ijn(1, 18, 3091.9688604755),
				ijn(1, 44, 252266.40357872),
				ijn(2, 0, -6.1707422868339E-03),
				ijn(2, 2, -0.31078046629583),
				ijn(2, 7, 11.670873077107),
				ijn(2, 36, 128127984.04046),
				ijn(2, 38, -985549096.23276),
				ijn(2, 40, 2822454697.3002),
				ijn(2, 42, -3594897141.0703),
				ijn(2, 44, 1722734991.3197),
				ijn(3, 24, -13551.334240775),
				ijn(3, 44, 12848734.66465),
				ijn(4, 12, 1.3865724283226),
				ijn(4, 32, 235988.32556514),
				ijn(4, 44, -13105236.545054),
				ijn(5, 32, 7399.9835474766),
				ijn(5, 36, -551966.9703006),
				ijn(5, 42, 3715408.5996233),
				ijn(6, 34, 19127.72923966),
				ijn(6, 44, -415351.64835634),
				ijn(7, 28, -62.459855192507)]
	tph_bw2b=[	ijn(0, 0, 1489.5041079516),
				ijn(0, 1, 743.07798314034),
				ijn(0, 2, -97.708318797837),
				ijn(0, 12, 2.4742464705674),
				ijn(0, 18, -0.63281320016026),
				ijn(0, 24, 1.1385952129658),
				ijn(0, 28, -0.47811863648625),
				ijn(0, 40, 8.5208123431544E-03),
				ijn(1, 0, 0.93747147377932),
				ijn(1, 2, 3.3593118604916),
				ijn(1, 6, 3.3809355601454),
				ijn(1, 12, 0.16844539671904),
				ijn(1, 18, 0.73875745236695),
				ijn(1, 24, -0.47128737436186),
				ijn(1, 28, 0.15020273139707),
				ijn(1, 40, -0.002176411421975),
				ijn(2, 2, -0.021810755324761),
				ijn(2, 8, -0.10829784403677),
				ijn(2, 18, -0.046333324635812),
				ijn(2, 40, 7.1280351959551E-05),
				ijn(3, 1, 1.1032831789999E-04),
				ijn(3, 2, 1.8955248387902E-04),
				ijn(3, 12, 3.0891541160537E-03),
				ijn(3, 24, 1.3555504554949E-03),
				ijn(4, 2, 2.8640237477456E-07),
				ijn(4, 12, -1.0779857357512E-05),
				ijn(4, 18, -7.6462712454814E-05),
				ijn(4, 24, 1.4052392818316E-05),
				ijn(4, 28, -3.1083814331434E-05),
				ijn(4, 40, -1.0302738212103E-06),
				ijn(5, 18, 2.821728163504E-07),
				ijn(5, 24, 1.2704902271945E-06),
				ijn(5, 40, 7.3803353468292E-08),
				ijn(6, 28, -1.1030139238909E-08),
				ijn(7, 2, -8.1456365207833E-14),
				ijn(7, 28, -2.5180545682962E-11),
				ijn(9, 1, -1.7565233969407E-18),
				ijn(9, 40, 8.6934156344163E-15)]						
	tph_bw2c=[	ijn(-7, 0, -3236839855524.2),
				ijn(-7, 4, 7326335090218.1),
				ijn(-6, 0, 358250899454.47),
				ijn(-6, 2, -583401318515.9),
				ijn(-5, 0, -10783068217.47),
				ijn(-5, 2, 20825544563.171),
				ijn(-2, 0, 610747.83564516),
				ijn(-2, 1, 859777.2253558),
				ijn(-1, 0, -25745.72360417),
				ijn(-1, 2, 31081.088422714),
				ijn(0, 0, 1208.2315865936),
				ijn(0, 1, 482.19755109255),
				ijn(1, 4, 3.7966001272486),
				ijn(1, 8, -10.842984880077),
				ijn(2, 4, -0.04536417267666),
				ijn(6, 0, 1.4559115658698E-13),
				ijn(6, 1, 1.126159740723E-12),
				ijn(6, 4, -1.7804982240686E-11),
				ijn(6, 10, 1.2324579690832E-07),
				ijn(6, 12, -1.1606921130984E-06),
				ijn(6, 16, 2.7846367088554E-05),
				ijn(6, 20, -5.9270038474176E-04),
				ijn(6, 22, 1.2918582991878E-03)]
	tps_bw2a = [fijn(-1.50, -24, -3.9235983861984E+05),
				fijn(-1.50, -23, 5.1526573827270E+05),
				fijn(-1.50, -19, 4.0482443161048E+04),
				fijn(-1.50, -13, -3.2193790923902E+02),
				fijn(-1.50, -11, 9.6961424218694E+01),
				fijn(-1.50, -10, -2.2867846371773E+01),
				fijn(-1.25, -19, -4.4942914124357E+05),
				fijn(-1.25, -15, -5.0118336020166E+03),
				fijn(-1.25, -6, 3.5684463560015E-01),
				fijn(-1.00, -26, 4.4235335848190E+04),
				fijn(-1.00, -21, -1.3673388811708E+04),
				fijn(-1.00, -17, 4.2163260207864E+05),
				fijn(-1.00, -16, 2.2516925837475E+04),
				fijn(-1.00, -9, 4.7442144865646E+02),
				fijn(-1.00, -8, -1.4931130797647E+02),
				fijn(-0.75, -15, -1.9781126320452E+05),
				fijn(-0.75, -14, -2.3554399470760E+04),
				fijn(-0.50, -26, -1.9070616302076E+04),
				fijn(-0.50, -13, 5.5375669883164E+04),
				fijn(-0.50, -9, 3.8293691437363E+03),
				fijn(-0.50, -7, -6.0391860580567E+02),
				fijn(-0.25, -27, 1.9363102620331E+03),
				fijn(-0.25, -25, 4.2660643698610E+03),
				fijn(-0.25, -11, -5.9780638872718E+03),
				fijn(-0.25, -6, -7.0401463926862E+02),
				fijn(0.25, 1, 3.3836784107553E+02),
				fijn(0.25, 4, 2.0862786635187E+01),
				fijn(0.25, 8, 3.3834172656196E-02),
				fijn(0.25, 11, -4.3124428414893E-05),
				fijn(0.50, 0, 1.6653791356412E+02),
				fijn(0.50, 1, -1.3986292055898E+02),
				fijn(0.50, 5, -7.8849547999872E-01),
				fijn(0.50, 6, 7.2132411753872E-02),
				fijn(0.50, 10, -5.9754839398283E-03),
				fijn(0.50, 14, -1.2141358953904E-05),
				fijn(0.50, 16, 2.3227096733871E-07),
				fijn(0.75, 0, -1.0538463566194E+01),
				fijn(0.75, 4, 2.0718925496502E+00),
				fijn(0.75, 9, -7.2193155260427E-02),
				fijn(0.75, 17, 2.0749887081120E-07),
				fijn(1.00, 7, -1.8340657911379E-02),
				fijn(1.00, 18, 2.9036272348696E-07),
				fijn(1.25, 3, 2.1037527893619E-01),
				fijn(1.25, 15, 2.5681239729999E-04),
				fijn(1.50, 5, -1.2799002933781E-02),
				fijn(1.50, 18, -8.2198102652018E-06)]
	tps_bw2b = [ijn(-6, 0, 3.1687665083497E+05),
				ijn(-6, 11, 2.0864175881858E+01),
				ijn(-5, 0, -3.9859399803599E+05),
				ijn(-5, 11, -2.1816058518877E+01),
				ijn(-4, 0, 2.2369785194242E+05),
				ijn(-4, 1, -2.7841703445817E+03),
				ijn(-4, 11, 9.9207436071480E+00),
				ijn(-3, 0, -7.5197512299157E+04),
				ijn(-3, 1, 2.9708605951158E+03),
				ijn(-3, 11, -3.4406878548526E+00),
				ijn(-3, 12, 3.8815564249115E-01),
				ijn(-2, 0, 1.7511295085750E+04),
				ijn(-2, 1, -1.4237112854449E+03),
				ijn(-2, 6, 1.0943803364167E+00),
				ijn(-2, 10, 8.9971619308495E-01),
				ijn(-1, 0, -3.3759740098958E+03),
				ijn(-1, 1, 4.7162885818355E+02),
				ijn(-1, 5, -1.9188241993679E+00),
				ijn(-1, 8, 4.1078580492196E-01),
				ijn(-1, 9, -3.3465378172097E-01),
				ijn(0, 0, 1.3870034777505E+03),
				ijn(0, 1, -4.0663326195838E+02),
				ijn(0, 2, 4.1727347159610E+01),
				ijn(0, 4, 2.1932549434532E+00),
				ijn(0, 5, -1.0320050009077E+00),
				ijn(0, 6, 3.5882943516703E-01),
				ijn(0, 9, 5.2511453726066E-03),
				ijn(1, 0, 1.2838916450705E+01),
				ijn(1, 1, -2.8642437219381E+00),
				ijn(1, 2, 5.6912683664855E-01),
				ijn(1, 3, -9.9962954584931E-02),
				ijn(1, 7, -3.2632037778459E-03),
				ijn(1, 8, 2.3320922576723E-04),
				ijn(2, 0, -1.5334809857450E-01),
				ijn(2, 1, 2.9072288239902E-02),
				ijn(2, 5, 3.7534702741167E-04),
				ijn(3, 0, 1.7296691702411E-03),
				ijn(3, 1, -3.8556050844504E-04),
				ijn(3, 3, -3.5017712292608E-05),
				ijn(4, 0, -1.4566393631492E-05),
				ijn(4, 1, 5.6420857267269E-06),
				ijn(5, 0, 4.1286150074605E-08),
				ijn(5, 1, -2.0684671118824E-08),
				ijn(5, 2, 1.6409393674725E-09)]
	tps_bw2c = [ijn(-2, 0, 9.0968501005365E+02),
				ijn(-2, 1, 2.4045667088420E+03),
				ijn(-1, 0, -5.9162326387130E+02),
				ijn(0, 0, 5.4145404128074E+02),
				ijn(0, 1, -2.7098308411192E+02),
				ijn(0, 2, 9.7976525097926E+02),
				ijn(0, 3, -4.6966772959435E+02),
				ijn(1, 0, 1.4399274604723E+01),
				ijn(1, 1, -1.9104204230429E+01),
				ijn(1, 3, 5.3299167111971E+00),
				ijn(1, 4, -2.1252975375934E+01),
				ijn(2, 0, -3.1147334413760E-01),
				ijn(2, 1, 6.0334840894623E-01),
				ijn(2, 2, -4.2764839702509E-02),
				ijn(3, 0, 5.8185597255259E-03),
				ijn(3, 1, -1.4597008284753E-02),
				ijn(3, 5, 5.6631175631027E-03),
				ijn(4, 0, -7.6155864584577E-05),
				ijn(4, 1, 2.2440342919332E-04),
				ijn(4, 4, -1.2561095013413E-05),
				ijn(5, 0, 6.3323132660934E-07),
				ijn(5, 1, -2.0541989675375E-06),
				ijn(5, 2, 3.6405370390082E-08),
				ijn(6, 0, -2.9759897789215E-09),
				ijn(6, 1, 1.0136618529763E-08),
				ijn(7, 0, 5.9925719692351E-12),
				ijn(7, 1, -2.0677870105164E-11),
				ijn(7, 3, -2.0874278181886E-11),
				ijn(7, 4, 1.0162166825089E-10),
				ijn(7, 5, -1.6429828281347E-10)]
	phs_bw2a = [ijn(0, 1, -1.82575361923032E-02),
				ijn(0, 3, -1.25229548799536E-01),
				ijn(0, 6, 5.92290437320145E-01),
				ijn(0, 16, 6.04769706185122E+00),
				ijn(0, 20, 2.38624965444474E+02),
				ijn(0, 22, -2.98639090222922E+02),
				ijn(1, 0, 5.12250813040750E-02),
				ijn(1, 1, -4.37266515606486E-01),
				ijn(1, 2, 4.13336902999504E-01),
				ijn(1, 3, -5.16468254574773E+00),
				ijn(1, 5, -5.57014838445711E+00),
				ijn(1, 6, 1.28555037824478E+01),
				ijn(1, 10, 1.14144108953290E+01),
				ijn(1, 16, -1.19504225652714E+02),
				ijn(1, 20, -2.84777985961560E+03),
				ijn(1, 22, 4.31757846408006E+03),
				ijn(2, 3, 1.12894040802650E+00),
				ijn(2, 16, 1.97409186206319E+03),
				ijn(2, 20, 1.51612444706087E+03),
				ijn(3, 0, 1.41324451421235E-02),
				ijn(3, 2, 5.85501282219601E-01),
				ijn(3, 3, -2.97258075863012E+00),
				ijn(3, 6, 5.94567314847319E+00),
				ijn(3, 16, -6.23656565798905E+03),
				ijn(4, 16, 9.65986235133332E+03),
				ijn(5, 3, 6.81500934948134E+00),
				ijn(5, 16, -6.33207286824489E+03),
				ijn(6, 3, -5.58919224465760E+00),
				ijn(7, 1, 4.00645798472063E-02)]
	phs_bw2b = [ijn(0, 0, 8.01496989929495E-02),
				ijn(0, 1, -5.43862807146111E-01),
				ijn(0, 2, 3.37455597421283E-01),
				ijn(0, 4, 8.90555451157450E+00),
				ijn(0, 8, 3.13840736431485E+02),
				ijn(1, 0, 7.97367065977789E-01),
				ijn(1, 1, -1.21616973556240E+00),
				ijn(1, 2, 8.72803386937477E+00),
				ijn(1, 3, -1.69769781757602E+01),
				ijn(1, 5, -1.86552827328416E+02),
				ijn(1, 12, 9.51159274344237E+04),
				ijn(2, 1, -1.89168510120494E+01),
				ijn(2, 6, -4.33407037194840E+03),
				ijn(2, 18, 5.43212633012715E+08),
				ijn(3, 0, 1.44793408386013E-01),
				ijn(3, 1, 1.28024559637516E+02),
				ijn(3, 7, -6.72309534071268E+04),
				ijn(3, 12, 3.36972380095287E+07),
				ijn(4, 1, -5.86634196762720E+02),
				ijn(4, 16, -2.21403224769889E+10),
				ijn(5, 1, 1.71606668708389E+03),
				ijn(5, 12, -5.70817595806302E+08),
				ijn(6, 1, -3.12109693178482E+03),
				ijn(6, 8, -2.07841384633010E+06),
				ijn(6, 18, 3.05605946157786E+12),
				ijn(7, 1, 3.22157004314333E+03),
				ijn(7, 16, 3.26810259797295E+11),
				ijn(8, 1, -1.44104158934487E+03),
				ijn(8, 3, 4.10694867802691E+02),
				ijn(8, 14, 1.09077066873024E+11),
				ijn(8, 18, -2.47964654258893E+13),
				ijn(12, 10, 1.88801906865134E+09),
				ijn(14, 16, -1.23651009018773E+14)]
	phs_bw2c = [ijn(0, 0, 1.12225607199012E-01),
				ijn(0, 1, -3.39005953606712E+00),
				ijn(0, 2, -3.20503911730094E+01),
				ijn(0, 3, -1.97597305104900E+02),
				ijn(0, 4, -4.07693861553446E+02),
				ijn(0, 8, 1.32943775222331E+04),
				ijn(1, 0, 1.70846839774007E+00),
				ijn(1, 2, 3.73694198142245E+01),
				ijn(1, 5, 3.58144365815434E+03),
				ijn(1, 8, 4.23014446424664E+05),
				ijn(1, 14, -7.51071025760063E+08),
				ijn(2, 2, 5.23446127607898E+01),
				ijn(2, 3, -2.28351290812417E+02),
				ijn(2, 7, -9.60652417056937E+05),
				ijn(2, 10, -8.07059292526074E+07),
				ijn(2, 18, 1.62698017225669E+12),
				ijn(3, 0, 7.72465073604171E-01),
				ijn(3, 5, 4.63929973837746E+04),
				ijn(3, 8, -1.37317885134128E+07),
				ijn(3, 16, 1.70470392630512E+12),
				ijn(3, 18, -2.51104628187308E+13),
				ijn(4, 18, 3.17748830835520E+13),
				ijn(5, 1, 5.38685623675312E+01),
				ijn(5, 4, -5.53089094625169E+04),
				ijn(5, 6, -1.02861522421405E+06),
				ijn(5, 14, 2.04249418756234E+12),
				ijn(6, 8, 2.73918446626977E+08),
				ijn(6, 18, -2.63963146312685E+15),
				ijn(10, 7, -1.07890854108088E+09),
				ijn(12, 7, -2.96492620980124E+10),
				ijn(16, 10, -1.11754907323424E+15)]
	tph_bw3a=[	ijn(-12, 0, -1.33645667811215E-07),
				ijn(-12, 1, 4.55912656802978E-06),
				ijn(-12, 2, -1.46294640700979E-05),
				ijn(-12, 6, 6.3934131297008E-03),
				ijn(-12, 14, 372.783927268847),
				ijn(-12, 16, -7186.54377460447),
				ijn(-12, 20, 573494.7521034),
				ijn(-12, 22, -2675693.29111439),
				ijn(-10, 1, -3.34066283302614E-05),
				ijn(-10, 5, -2.45479214069597E-02),
				ijn(-10, 12, 47.8087847764996),
				ijn(-8, 0, 7.64664131818904E-06),
				ijn(-8, 2, 1.28350627676972E-03),
				ijn(-8, 4, 1.71219081377331E-02),
				ijn(-8, 10, -8.51007304583213),
				ijn(-5, 2, -1.36513461629781E-02),
				ijn(-3, 0, -3.84460997596657E-06),
				ijn(-2, 1, 3.37423807911655E-03),
				ijn(-2, 3, -0.551624873066791),
				ijn(-2, 4, 0.72920227710747),
				ijn(-1, 0, -9.92522757376041E-03),
				ijn(-1, 2, -0.119308831407288),
				ijn(0, 0, 0.793929190615421),
				ijn(0, 1, 0.454270731799386),
				ijn(1, 1, 0.20999859125991),
				ijn(3, 0, -6.42109823904738E-03),
				ijn(3, 1, -0.023515586860454),
				ijn(4, 0, 2.52233108341612E-03),
				ijn(4, 3, -7.64885133368119E-03),
				ijn(10, 4, 1.36176427574291E-02),
				ijn(12, 5, -1.33027883575669E-02)]
	tph_bw3b=[	ijn(-12, 0, 3.2325457364492E-05),
				ijn(-12, 1, -1.27575556587181E-04),
				ijn(-10, 0, -4.75851877356068E-04),
				ijn(-10, 1, 1.56183014181602E-03),
				ijn(-10, 5, 0.105724860113781),
				ijn(-10, 10, -85.8514221132534),
				ijn(-10, 12, 724.140095480911),
				ijn(-8, 0, 2.96475810273257E-03),
				ijn(-8, 1, -5.92721983365988E-03),
				ijn(-8, 2, -1.26305422818666E-02),
				ijn(-8, 4, -0.115716196364853),
				ijn(-8, 10, 84.9000969739595),
				ijn(-6, 0, -1.08602260086615E-02),
				ijn(-6, 1, 1.54304475328851E-02),
				ijn(-6, 2, 7.50455441524466E-02),
				ijn(-4, 0, 2.52520973612982E-02),
				ijn(-4, 1, -6.02507901232996E-02),
				ijn(-3, 5, -3.07622221350501),
				ijn(-2, 0, -5.74011959864879E-02),
				ijn(-2, 4, 5.03471360939849),
				ijn(-1, 2, -0.925081888584834),
				ijn(-1, 4, 3.91733882917546),
				ijn(-1, 6, -77.314600713019),
				ijn(-1, 10, 9493.08762098587),
				ijn(-1, 14, -1410437.19679409),
				ijn(-1, 16, 8491662.30819026),
				ijn(0, 0, 0.861095729446704),
				ijn(0, 2, 0.32334644281172),
				ijn(1, 1, 0.873281936020439),
				ijn(3, 1, -0.436653048526683),
				ijn(5, 1, 0.286596714529479),
				ijn(6, 1, -0.131778331276228),
				ijn(8, 1, 6.76682064330275E-03)]
	vph_bw3a=[	ijn(-12, 6, 5.29944062966028E-03),
				ijn(-12, 8, -0.170099690234461),
				ijn(-12, 12, 11.1323814312927),
				ijn(-12, 18, -2178.98123145125),
				ijn(-10, 4, -5.06061827980875E-04),
				ijn(-10, 7, 0.556495239685324),
				ijn(-10, 10, -9.43672726094016),
				ijn(-8, 5, -0.297856807561527),
				ijn(-8, 12, 93.9353943717186),
				ijn(-6, 3, 1.92944939465981E-02),
				ijn(-6, 4, 0.421740664704763),
				ijn(-6, 22, -3689141.2628233),
				ijn(-4, 2, -7.37566847600639E-03),
				ijn(-4, 3, -0.354753242424366),
				ijn(-3, 7, -1.99768169338727),
				ijn(-2, 3, 1.15456297059049),
				ijn(-2, 16, 5683.6687581596),
				ijn(-1, 0, 8.08169540124668E-03),
				ijn(-1, 1, 0.172416341519307),
				ijn(-1, 2, 1.04270175292927),
				ijn(-1, 3, -0.297691372792847),
				ijn(0, 0, 0.560394465163593),
				ijn(0, 1, 0.275234661176914),
				ijn(1, 0, -0.148347894866012),
				ijn(1, 1, -6.51142513478515E-02),
				ijn(1, 2, -2.92468715386302),
				ijn(2, 0, 6.64876096952665E-02),
				ijn(2, 2, 3.52335014263844),
				ijn(3, 0, -1.46340792313332E-02),
				ijn(4, 2, -2.24503486668184),
				ijn(5, 2, 1.10533464706142),
				ijn(8, 2, -4.08757344495612E-02)]
	phs_bw3a = [ijn(0, 0, 7.70889828326934E+00),
				ijn(0, 1, -2.60835009128688E+01),
				ijn(0, 5, 2.67416218930389E+02),
				ijn(1, 0, 1.72221089496844E+01),
				ijn(1, 3, -2.93542332145970E+02),
				ijn(1, 4, 6.14135601882478E+02),
				ijn(1, 8, -6.10562757725674E+04),
				ijn(1, 14, -6.51272251118219E+07),
				ijn(2, 6, 7.35919313521937E+04),
				ijn(2, 16, -1.16646505914191E+10),
				ijn(3, 0, 3.55267086434461E+01),
				ijn(3, 2, -5.96144543825955E+02),
				ijn(3, 3, -4.75842430145708E+02),
				ijn(4, 0, 6.96781965359503E+01),
				ijn(4, 1, 3.35674250377312E+02),
				ijn(4, 4, 2.50526809130882E+04),
				ijn(4, 5, 1.46997380630766E+05),
				ijn(5, 28, 5.38069315091534E+19),
				ijn(6, 28, 1.43619827291346E+21),
				ijn(7, 24, 3.64985866165994E+19),
				ijn(8, 1, -2.54741561156775E+03),
				ijn(10, 32, 2.40120197096563E+27),
				ijn(10, 36, -3.93847464679496E+29),
				ijn(14, 22, 1.47073407024852E+24),
				ijn(18, 28, -4.26391250432059E+31),
				ijn(20, 36, 1.94509340621077E+38),
				ijn(22, 16, 6.66212132114896E+23),
				ijn(22, 28, 7.06777016552858E+33),
				ijn(24, 36, 1.75563621975576E+41),
				ijn(28, 16, 1.08408607429124E+28),
				ijn(28, 36, 7.30872705175151E+43),
				ijn(32, 10, 1.59145847398870E+24),
				ijn(32, 28, 3.77121605943324E+40)]
	phs_bw3b = [ijn(-12, 2, 1.25244360717979E-13),
				ijn(-12, 10, -1.26599322553713E-02),
				ijn(-12, 12, 5.06878030140626E+00),
				ijn(-12, 14, 3.17847171154202E+01),
				ijn(-12, 20, -3.91041161399932E+05),
				ijn(-10, 2, -9.75733406392044E-11),
				ijn(-10, 10, -1.86312419488279E+01),
				ijn(-10, 14, 5.10973543414101E+02),
				ijn(-10, 18, 3.73847005822362E+05),
				ijn(-8, 2, 2.99804024666572E-08),
				ijn(-8, 8, 2.00544393820342E+01),
				ijn(-6, 2, -4.98030487662829E-06),
				ijn(-6, 6, -1.02301806360030E+01),
				ijn(-6, 7, 5.52819126990325E+01),
				ijn(-6, 8, -2.06211367510878E+02),
				ijn(-5, 10, -7.94012232324823E+03),
				ijn(-4, 4, 7.82248472028153E+00),
				ijn(-4, 5, -5.86544326902468E+01),
				ijn(-4, 8, 3.55073647696481E+03),
				ijn(-3, 1, -1.15303107290162E-04),
				ijn(-3, 3, -1.75092403171802E+00),
				ijn(-3, 5, 2.57981687748160E+02),
				ijn(-3, 6, -7.27048374179467E+02),
				ijn(-2, 0, 1.21644822609198E-04),
				ijn(-2, 1, 3.93137871762692E-02),
				ijn(-1, 0, 7.04181005909296E-03),
				ijn(0, 3, -8.29108200698110E+01),
				ijn(2, 0, -2.65178818131250E-01),
				ijn(2, 1, 1.37531682453991E+01),
				ijn(5, 0, -5.22394090753046E+01),
				ijn(6, 1, 2.40556298941048E+03),
				ijn(8, 1, -2.27361631268929E+04),
				ijn(10, 1, 8.90746343932567E+04),
				ijn(14, 3, -2.39234565822486E+07),
				ijn(14, 7, 5.68795808129714E+09)]	
	ths_bw4 = [	ijn(0, 0, 1.79882673606601E-01),
				ijn(0, 3, -2.67507455199603E-01),
				ijn(0, 12, 1.16276722612600E+00),
				ijn(1, 0, 1.47545428713616E-01),
				ijn(1, 1, -5.12871635973248E-01),
				ijn(1, 2, 4.21333567697984E-01),
				ijn(1, 5, 5.63749522189870E-01),
				ijn(2, 0, 4.29274443819153E-01),
				ijn(2, 5, -3.35704552142140E+00),
				ijn(2, 8, 1.08890916499278E+01),
				ijn(3, 0, -2.48483390456012E-01),
				ijn(3, 2, 3.04153221906390E-01),
				ijn(3, 3, -4.94819763939905E-01),
				ijn(3, 4, 1.07551674933261E+00),
				ijn(4, 0, 7.33888415457688E-02),
				ijn(4, 1, 1.40170545411085E-02),
				ijn(5, 1, -1.06110975998808E-01),
				ijn(5, 2, 1.68324361811875E-02),
				ijn(5, 4, 1.25028363714877E+00),
				ijn(5, 16, 1.01316840309509E+03),
				ijn(6, 6, -1.51791558000712E+00),
				ijn(6, 8, 5.24277865990866E+01),
				ijn(6, 22, 2.30495545563912E+04),
				ijn(8, 1, 2.49459806365456E-02),
				ijn(10, 20, 2.10796467412137E+06),
				ijn(10, 36, 3.66836848613065E+08),
				ijn(12, 24, -1.44814105365163E+08),
				ijn(14, 1, -1.79276373003590E-03),
				ijn(14, 28, 4.89955602100459E+09),
				ijn(16, 12, 4.71262212070518E+02),
				ijn(16, 32, -8.29294390198652E+10),
				ijn(18, 14, -1.71545662263191E+03),
				ijn(18, 22, 3.55777682973575E+06),
				ijn(18, 36, 5.86062760258436E+11),
				ijn(20, 24, -1.29887635078195E+07),
				ijn(28, 36, 3.17247449371057E+10)]
)


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

// P=f(h,s), region 3
fn bw3_p_hs(h f64, s f64) f64{
    sc := 4.41202148223476
    if s <= sc{
        return bw3a_p_hs(h, s)
	}else{
        return bw3b_p_hs(h, s)
	}
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
  

