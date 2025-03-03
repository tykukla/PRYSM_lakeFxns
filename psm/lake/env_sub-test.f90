      program lakemodel

      implicit none
      ! Version 1.0
! Declarations ********************************************************
      real dz,dt,readdt,basedep,zo,z_screen,delta,rhowat,rhoice,rhosnow
      real Le,Lei,fusion,surf,fracmin,fraclim,qwtau,cpw_ice,condw,beta
      real snocrit,kv,dm,pi,tsurfp,mixing,trace,depth_begin,salty_begin
      real trace_i,trace_a,d_frac,tempi_a,hice_a,hsnow_a,xlat,xlon,gmt
      real b_area,runout_sum,alb_slush,alb_snow,fsol1,fsol2,sical0
      real sical1
      real salty_a,eta,fraci_a,temp_a,d_fraca,snowcut,o18air,f,ti_a
      real deutair,dexch,alphak,lamisw,lamilw,ca,cb,c70,c71,c72,c73
      real lamssw,lamslw,afrac1,afrac2,condi,conds,cpair,emis
      real o18snow_a,deutsnow_a,psum_a,surfarea_a,hsnow_ave,temp_ave
      real tsurf_ave,fice_ave,evap_ave,hice_ave,o18_ave,deut_ave
      real oblq,cdrn,raddeg,dpd,grav,sigma,rair,rvap,cvap
      real co18prec,co18run,cdeutprec,cdeutrun,xt,mix_ave

      integer max_dep,ix1,iy1,n_trace,i_area,lcount,nspin
      integer iin,iout,numpts,num_tra,ktau,iupto,ilake,jlake
      integer depth,depth_a,ktauwan,mixmax_a

      logical bndry_flag,iceflag,s_flag,o18flag,deutflag,ar_flag
      logical snow_flag_a,wb_flag,melt_flag_a

      real tempinit, deutinit, o18init !Ashling
!**********************************************************************
! Lake specific parameters ********************************************

      real, parameter :: oblq = 23.4           ! obliquity
      real, parameter :: xlat = -6.30          ! latitude (negative for South)
      real, parameter :: xlon = 29.5           ! longitude (negative for West)
      real, parameter :: gmt = +3              ! local time relative to gmt in hours
      real, parameter :: max_dep = 999         ! depth of lake at sill in meters 
      real, parameter :: basedep = 733.	    ! elevation of basin bottom in meters
      real, parameter :: b_area = 23100000.    ! area of catchment+lake in hectares 
      real, parameter :: cdrn = 2.0e-3         ! neutral drag coefficient 1.8 HAD 1.7GISS 1.2CCSM
      real, parameter :: eta = 0.065             ! shortwave extinction coefficient (1/m)
      real, parameter :: f= 0.3                 ! fraction of advected air over lake
      real, parameter :: alb_slush = 0.4       ! albedo of melting snow
      real, parameter :: alb_snow = 0.7        ! albedo of non-melting snow

      real, parameter :: depth_begin = 570       ! prescribed depth in meters
      real, parameter :: salty_begin = 0.0     ! prescribed salinity in ppt

      real, parameter :: o18air = -14.0         ! d18O of air above lake
      real, parameter :: deutair = -96.       ! dD of air above lake


      ! Lake area moved to bottom

      !real area(max_dep)
      !data area/570*3290000./           ! lake area in hectares by depth
     
      !integer, parameter :: hyplength=1500
      !real, dimension(hyplength):: area
      !open(unit=100, file='hypsograph.txt') 
      !read(100,*) area
      !close(unit=100) ! for orderly closing

      !Begin Ashling
      ! Removed from the model, for easier access
      character(38) :: datafile='ERA_INTERIM_climatology_Tanganyika.txt' ! the data file to open in FILE_OPEN subroutine
      real, parameter :: tempinit= 23.0 ! temperature to initialize lake at in INIT_LAKE subroutine
      real, parameter :: deutinit= 24.0! dD to initialize lake at in INIT_LAKE subroutine
      real, parameter :: o18init= 3.7  ! d18O to initialize lake at in INIT_LAKE subroutine
      !End Ashling

! Simulation specific parameters **************************************

      real, parameter :: nspin = 30              ! number of years for spinup
      real, parameter :: bndry_flag = .false.   ! true for explict boundry layer computations; presently only for sigma coord climate models
      real, parameter :: sigma = 0.9925561       ! sigma level for boundary flag
      real, parameter :: wb_flag = .true.      ! true for variable lake depth
      real, parameter :: iceflag = .false.        ! true for variable ice cover
      real, parameter :: s_flag = .false.        ! true for variable salinity
      real, parameter :: o18flag = .true.       ! true for variable d18O 
      real, parameter :: deutflag = .true.      ! true for variable dD
      real, parameter :: z_screen = 5.0          ! height of met inputs

!**********************************************************************
! Other parameters DO NOT CHANGE without good reason for doing so******

      real, parameter :: dz=1.0                ! vertical layer thickness in m
      real, parameter :: dt = 1.*60.*30.       ! model time step in seconds
      real, parameter :: zo=0.0004             ! water surface roughness length
      real, parameter :: n_trace=3             ! model does 3 tracers (o18, D, sal)
      real, parameter :: delta=5.67e-8         ! s-b constant
      real, parameter :: rhowat = 1000.        ! density of water
      real, parameter :: rhosnow = 330.        ! density of snow
      real, parameter :: rhoice = 917.         ! density of ice
      real, parameter :: snowcut = 2.2         ! temp at which precip is snow
      real, parameter :: cpair = 1.00464e3     ! specific heat capacity dry air
      real, parameter :: rair = 287.05         ! specific gas constant dry air
      real, parameter :: rvap = 461.495        ! specific gas constant water vapor
      real, parameter :: cvap = 1.81e3         ! specific heat capacity water vapor
      real, parameter :: Le = 2.45e6           ! latent heat of vaporization water
      real, parameter :: Lei = 2.5e6           ! latent heat of vaporization ice
      real, parameter :: ca = 6.11e2           ! constant for Teten's formula
      real, parameter :: cb = 273.16           ! constant for Teten's formula
      real, parameter :: c72 = 17.269          ! constant for Teten over water
      real, parameter :: c73 = 35.86           ! constant for Teten over water
      real, parameter :: c70 = 21.874          ! constant for Teten over ice
      real, parameter :: c71 = 7.66            ! constant for Teten over ice
      real, parameter :: emis = 0.97           ! longwave emmisivity
      real, parameter :: fusion=3.34e5         ! latent heat of fusion
      real, parameter :: surf = 1.0            ! surface thickness over which the flux exchange occurs in water
      real, parameter :: fracmin= 0.01         ! min ice thick in meters
      real, parameter :: fraclim = 0.02        ! min ice fraction
      real, parameter :: qwtau = 86400.        ! D. Pollard sub-ice time constant
      real, parameter :: cpw_ice = 4200.       ! heat capacity of ice
      real, parameter :: condw = .58           ! conduction coefficinet for water
      real, parameter :: beta=0.4              ! frac of solar rad abs in water surface layer
      real, parameter :: fsol1=0.5             ! sw weighting for surface albedo
      real, parameter :: fsol2=0.5             ! lw weighting for surface albedo
      real, parameter :: sical0=0.6            ! sw albedo for sea ice
      real, parameter :: sical1=0.4            ! lw albedo for sea ice
      real, parameter :: kv=0.4                ! vonkarman constant
      real, parameter :: dm=1.38889E-07        ! thermal molecular diffusivity of water
      real, parameter :: pi=3.141592654        !
      real, parameter :: snocrit = 0.05        ! for albedo, in m
      real, parameter :: dexch  =1.0           ! fixed depth (m) over which isotope exchange occurs
      real, parameter :: alphak = 0.994        ! fractionation coeff
      real, parameter :: lamisw = 3.0          ! extinction coef for sw thru ice
      real, parameter :: lamilw = 20.          ! extinction coef for lw thru ice
      real, parameter :: lamssw = 3.0          ! extinction coef for sw thru snow
      real, parameter :: lamslw = 20.          ! extinction coef for lw thru snow
      real, parameter :: afrac1 = 0.7          ! fraction of light in visible band for ice model
      real, parameter :: afrac2 = 0.3          ! fraction of light in infrared band for ice model
      real, parameter :: condi = 2.3           ! thermal conductivity of ice
      real, parameter :: conds = 0.31          ! thermal conductivity of snow
      real, parameter :: raddeg = 0.0174533    ! radians per degree
      real, parameter :: dpd = 360./365.       ! degrees per day
      real, parameter :: grav = 9.80616        ! gravity

!**********************************************************************

! Common definitions
      common /ldata/ lcount, iin, iout, numpts , num_tra, ktau, iupto
      common mixing, o18snow_a,deutsnow_a
      common psum_a, depth_a
      common trace(max_dep,n_trace), trace_i(max_dep,n_trace)
      common trace_a (max_dep, n_trace),surfarea_a
      common d_fraca, tempi_a, hice_a, hsnow_a
      common salty_a, mixmax_a, fraci_a 
      common temp_a (max_dep), snow_flag_a, melt_flag_a
      common ktauwan (12),ti_a(max_dep)
      common mix_ave, tsurf_ave, fice_ave
      common evap_ave, hice_ave, hsnow_ave, o18_ave
      common deut_ave, runout_sum, temp_ave(max_dep)

!**********************************************************************
      !Ashling
      real, dimension(999) :: area
      data area/3290000,3246809,3235562,  &
      3227075,3219947,3213619,  &
      3207911,3202743,3197948,  &
      3193338,3188901,3184763,  &
      3180759,3176884,3173213,  &
      3169402,3165770,3162305,  &
      3158823,3155504,3151633,  &
      3147965,3144780,3141458,  &
      3138169,3134920,3131680,  &
      3128600,3125470,3122458,  &
      3119361,3116264,3113250,  &
      3110231,3107155,3104048,  &
      3100808,3097635,3094222,  &
      3090791,3086891,3083211,  &
      3079875,3076551,3073218,  &
      3069895,3066548,3063255,  &
      3059773,3056249,3051980,  &
      3047720,3044308,3041014,  &
      3037787,3034670,3031583,  &
      3028538,3025529,3022673,  &
      3019459,3016373,3013529,  &
      3010694,3007868,3005091,  &
      3002249,2999488,2996796,  &
      2993993,2991234,2988320,  &
      2985593,2982797,2980008,  &
      2977115,2974295,2971594,  &
      2968793,2966023,2962862,  &
      2959644,2956914,2954219,  &
      2951607,2948952,2946271,  &
      2943708,2941094,2938531,  &
      2935822,2933124,2930490,  &
      2927827,2925134,2922339,  &
      2919593,2916844,2914098,  &
      2911225,2907195,2903030,  &
      2900256,2897525,2894890,  &
      2892295,2889543,2886990,  &
      2884535,2882136,2879505,  &
      2876973,2874619,2872208,  &
      2869853,2867386,2864927,  &
      2862476,2860053,2857578,  &
      2854812,2852016,2849537,  &
      2847005,2844493,2841816,  &
      2839093,2836582,2834102,  &
      2831613,2828976,2826403,  &
      2823847,2821294,2818706,  &
      2815978,2813257,2810534,  &
      2807774,2804988,2802002,  &
      2798948,2795989,2792856,  &
      2789736,2786290,2782692,  &
      2779057,2775216,2771149,  &
      2765965,2761472,2758121,  &
      2754814,2751678,2748429,  &
      2745299,2742340,2739416,  &
      2736525,2733312,2730198,  &
      2727338,2724493,2721624,  &
      2718664,2715843,2713137,  &
      2710451,2707827,2705074,  &
      2702320,2699609,2697023,  &
      2694307,2691420,2688536,  &
      2685678,2682873,2679982,  &
      2676746,2673719,2670837,  &
      2667905,2664960,2662007,  &
      2659033,2656164,2653176,  &
      2650270,2647133,2643974,  &
      2640922,2637719,2634337,  &
      2630823,2627209,2623465,  &
      2619478,2614907,2608009,  &
      2602009,2597819,2593879,  &
      2590066,2586545,2583280,  &
      2580273,2577382,2574587,  &
      2571832,2569042,2566377,  &
      2563864,2561260,2558705,  &
      2556098,2553637,2551249,  &
      2548788,2546291,2543735,  &
      2541396,2539040,2536659,  &
      2534186,2531760,2529339,  &
      2526947,2524612,2522182,  &
      2519756,2517388,2514999,  &
      2512674,2510209,2507727,  &
      2505308,2502828,2500382,  &
      2497649,2494883,2492313,  &
      2489665,2486958,2484041,  &
      2481093,2478105,2474962,  &
      2471819,2467415,2462943,  &
      2459732,2456695,2453592,  &
      2450498,2447700,2445039,  &
      2442352,2439824,2437172,  &
      2434522,2432023,2429540,  &
      2427132,2424632,2422221,  &
      2419861,2417472,2415107,  &
      2412745,2410443,2408186,  &
      2405846,2403590,2401177,  &
      2398846,2396506,2394200,  &
      2391857,2389390,2387018,  &
      2384659,2382360,2379964,  &
      2377405,2374908,2372521,  &
      2370143,2367645,2365028,  &
      2362451,2359744,2357081,  &
      2354298,2351306,2348188,  &
      2344983,2341586,2337713,  &
      2330287,2321769,2317117,  &
      2313063,2309294,2305612,  &
      2302077,2298786,2295646,  &
      2292527,2289409,2286354,  &
      2283296,2280271,2277237,  &
      2274005,2270728,2267382,  &
      2264142,2261082,2258157,  &
      2255125,2252243,2249313,  &
      2246576,2243680,2240812,  &
      2238121,2235423,2232700,  &
      2230025,2227365,2224766,  &
      2222130,2219581,2217080,  &
      2214638,2212179,2209764,  &
      2207197,2204543,2201823,  &
      2199108,2196277,2193484,  &
      2190538,2187464,2184000,  &
      2180302,2176056,2169272,  &
      2162757,2158476,2154595,  &
      2150874,2147176,2143568,  &
      2140203,2136938,2133868,  &
      2130647,2127582,2124689,  &
      2121926,2119367,2116644,  &
      2113962,2111353,2108881,  &
      2106290,2103725,2101225,  &
      2098763,2096362,2093953,  &
      2091552,2089064,2086579,  &
      2084050,2081451,2078706,  &
      2075998,2073360,2070728,  &
      2068182,2065508,2062881,  &
      2060209,2057522,2054805,  &
      2052041,2049191,2046261,  &
      2043315,2040284,2037282,  &
      2034073,2030728,2027326,  &
      2023610,2018527,2013427,  &
      2009730,2006196,2002742,  &
      1999303,1996171,1993094,  &
      1990087,1987022,1983867,  &
      1980823,1977813,1974825,  &
      1971978,1969156,1966311,  &
      1963510,1960692,1957944,  &
      1955109,1952284,1949778,  &
      1947248,1944839,1942361,  &
      1939985,1937652,1935338,  &
      1932974,1930638,1928272,  &
      1925972,1923600,1921180,  &
      1918811,1916356,1913810,  &
      1911289,1908772,1906012,  &
      1903179,1900462,1897769,  &
      1895019,1892302,1889430,  &
      1886375,1883332,1880088,  &
      1875679,1871143,1867741,  &
      1864556,1861342,1858223,  &
      1855249,1852260,1849273,  &
      1846133,1842547,1838867,  &
      1835994,1833270,1830820,  &
      1828348,1825904,1823459,  &
      1821054,1818780,1816459,  &
      1814141,1811771,1809438,  &
      1807055,1804692,1802311,  &
      1799838,1797262,1794510,  &
      1791498,1788159,1785255,  &
      1782536,1779865,1777296,  &
      1774691,1772161,1769620,  &
      1767073,1764480,1761892,  &
      1759190,1756542,1753789,  &
      1750978,1747992,1744926,  &
      1741685,1738000,1732322,  &
      1727192,1723322,1719961,  &
      1716686,1713521,1710455,  &
      1707491,1704534,1701613,  &
      1698803,1695987,1693277,  &
      1690485,1687755,1684959,  &
      1682216,1679383,1676513,  &
      1673555,1670360,1667182,  &
      1664126,1661105,1658078,  &
      1655084,1652092,1649005,  &
      1645992,1642821,1639675,  &
      1636485,1633204,1629816,  &
      1626396,1623036,1619621,  &
      1616134,1612533,1608917,  &
      1604880,1600935,1596931,  &
      1592829,1588553,1584233,  &
      1579895,1575396,1570579,  &
      1565211,1556220,1547382,  &
      1540840,1534495,1528629,  &
      1523216,1518725,1514838,  &
      1511039,1507454,1503784,  &
      1500291,1497112,1494025,  &
      1491136,1488261,1485453,  &
      1482732,1480246,1477875,  &
      1475436,1473059,1470752,  &
      1468367,1466057,1463819,  &
      1461495,1459214,1456920,  &
      1454553,1452052,1449372,  &
      1447052,1444756,1442573,  &
      1440399,1438328,1436213,  &
      1434106,1432025,1429959,  &
      1427870,1425828,1423680,  &
      1421504,1419213,1416990,  &
      1414525,1411973,1409286,  &
      1405701,1402026,1399330,  &
      1396674,1394128,1391710,  &
      1389321,1387016,1384909,  &
      1382736,1380756,1378732,  &
      1376685,1374670,1372611,  &
      1370319,1367976,1365653,  &
      1363400,1360956,1358690,  &
      1356537,1354582,1352729,  &
      1350908,1349196,1347435,  &
      1345722,1344017,1342356,  &
      1340690,1339010,1337407,  &
      1335768,1334122,1332419,  &
      1330795,1329167,1327439,  &
      1325709,1323938,1322109,  &
      1320226,1318405,1316708,  &
      1314957,1313121,1311271,  &
      1309473,1307559,1305320,  &
      1302913,1301081,1299394,  &
      1297723,1295993,1294443,  &
      1292944,1291351,1289832,  &
      1288350,1286873,1285401,  &
      1283926,1282460,1281026,  &
      1279556,1278106,1276721,  &
      1275290,1273846,1272453,  &
      1271046,1269600,1268158,  &
      1266767,1265394,1263954,  &
      1262571,1261186,1259761,  &
      1258293,1256882,1255502,  &
      1254032,1252587,1251118,  &
      1249723,1248160,1246703,  &
      1245245,1243634,1242149,  &
      1240558,1238938,1237252,  &
      1235666,1233875,1232051,  &
      1230292,1227605,1225165,  &
      1223475,1221674,1219946,  &
      1218282,1216574,1214924,  &
      1213363,1211744,1210206,  &
      1208671,1207135,1205619,  &
      1204101,1202558,1201113,  &
      1199601,1198057,1196565,  &
      1195066,1193563,1192002,  &
      1190502,1188989,1187440,  &
      1185916,1184359,1182805,  &
      1181254,1179607,1178026,  &
      1176438,1174812,1173201,  &
      1171446,1169852,1168157,  &
      1166459,1164726,1163008,  &
      1161240,1159406,1157471,  &
      1155558,1153504,1151404,  &
      1149334,1146981,1144677,  &
      1141569,1138529,1136338,  &
      1133981,1131991,1129880,  &
      1127858,1125903,1123924,  &
      1122039,1120178,1118346,  &
      1116492,1114685,1112833,  &
      1111024,1109151,1107362,  &
      1105552,1103738,1101913,  &
      1100074,1098272,1096426,  &
      1094506,1092667,1090739,  &
      1088939,1087028,1085158,  &
      1083308,1081432,1079514,  &
      1077559,1075657,1073721,  &
      1071712,1069745,1067788,  &
      1065680,1063721,1061625,  &
      1059492,1057416,1055217,  &
      1053018,1050756,1048438,  &
      1046019,1043587,1040492,  &
      1037133,1034536,1031982,  &
      1029450,1026945,1024522,  &
      1022102,1019697,1017335,  &
      1014985,1012641,1010253,  &
      1007892,1005412,1002815,  &
      1000423,998109,995884,  &
      993705,991509,989445,  &
      987359,985307,983307,  &
      981316,979352,977579,  &
      976091,974587,973148,  &
      971682,970201,968740,  &
      967319,965858,964431,  &
      963090,961669,960262,  &
      958879,957493,956110,  &
      954653,953277,951890,  &
      950704,949531,948368,  &
      947163,945836,944400,  &
      943123,941817,940541,  &
      939325,938090,936929,  &
      935864,934755,933698,  &
      932679,931644,930610,  &
      929605,928601,927590,  &
      926607,925570,924599,  &
      923593,922584,921614,  &
      920593,919545,918549,  &
      917552,916502,915502,  &
      914476,913444,912403,  &
      911326,910224,909131,  &
      908073,907004,905888,  &
      904697,903544,902394,  &
      901124,899696,898292,  &
      896967,895786,894433,  &
      892921,891473,889954,  &
      887203,884429,882897,  &
      881498,880125,878838,  &
      877694,876537,875325,  &
      874286,873297,872265,  &
      871315,870364,869390,  &
      868426,867474,866485,  &
      865575,864630,863689,  &
      862718,861815,860843,  &
      859932,858968,857957,  &
      856985,856040,855072,  &
      854128,853098,852103,  &
      851073,850028,848962,  &
      847915,846838,845740,  &
      844585,843485,842302,  &
      841117,839878,838655,  &
      837315,835908,834481,  &
      832885,831156,828871,  &
      826733,825110,823530,  &
      822064,820675,819316,  &
      818038,816738,815477,  &
      814251,812961,811728,  &
      810515,809294,808126,  &
      806987,805768,804608,  &
      803416,802302,801183,  &
      800095,798953,797830,  &
      796764,795613,794478,  &
      793375,792284,791213,  &
      790079,788999,787901,  &
      786854,785689,784578,  &
      783426,782259,781059,  &
      779902,778748,777497,  &
      776261,775027,773730,  &
      772186,770693,769116/!,  &
      !767464/
      

      real year,day,ta_in,dp_in,ua_in,rlwd_in,sw_in,qa_in,  &
           prec_in,ps_in,runin_in,declin,ta_i,qa_i,ua_i,  &
           rh_i,sw_i,rlwd_i,ps_i,prec_i,runin_i,xtime,rh_in, &
           deutprec_in, deutrun_in, o18prec_in, o18run_in, &
           deutprec_i, deutrun_i, o18prec_i, o18run_i   !Ashling
      dimension year(2),day(2),ta_in(2),dp_in(2),ua_in(2),  &
                sw_in(2),rlwd_in(2),qa_in(2),prec_in(2),runin_in(2), &
                ps_in(2),rh_in(2), deutprec_in(2), deutrun_in(2), &
                o18prec_in(2), o18run_in(2)  !Ashling
      integer j,nsteps,ispin

      call file_open ! open input and output files
      call init_lake ! initialize lake variables
      ispin = 0

!	  Begin Ashling
      !Read in data based on 
      if (wb_flag) then
        !variable lake depth, with or without isotopes
        if(deutflag.and..not.o18flag) then
          !read for deut only
          read(15,*,end=998) year(1),day(1),ta_in(1),dp_in(1),   &
                         ua_in(1),sw_in(1),rlwd_in(1),ps_in(1),       &
                         prec_in(1), deutprec_in(1), &
                         runin_in(1), deutrun_in(1)
        else if(o18flag.and..not.deutflag) then
          !read for o18 only
          read(15,*,end=998) year(1),day(1),ta_in(1),dp_in(1),   &
                         ua_in(1),sw_in(1),rlwd_in(1),ps_in(1),       &
                         prec_in(1), o18prec_in(1), &
                         runin_in(1), o18run_in(1)          
        else if(o18flag.and.deutflag) then
          !read for deut and o18      
          read(15,*,end=998) year(1),day(1),ta_in(1),dp_in(1),   &
                         ua_in(1),sw_in(1),rlwd_in(1),ps_in(1),       &
                         prec_in(1), deutprec_in(1), o18prec_in(1), &
                         runin_in(1), deutrun_in(1), o18run_in(1)
        else
          !only wb_flag
           read(15,*,end=998) year(1),day(1),ta_in(1),dp_in(1),   &
                         ua_in(1),sw_in(1),rlwd_in(1),ps_in(1),       &
                         prec_in(1),runin_in(1)
        end if
      else
        !not variable lake depth, no isotopes
        read(15,*,end=998) year(1),day(1),ta_in(1),dp_in(1),   &
                         ua_in(1),sw_in(1),rlwd_in(1),ps_in(1)
        runin_in(1)=0.0
        prec_in(1)=0.0
      end if


      call datain (ta_in(1),dp_in(1),ua_in(1),sw_in(1),rlwd_in(1),  &
                   ps_in(1),prec_in(1),runin_in(1),qa_in(1),rh_in(1))

 150  if (wb_flag) then
        !variable lake depth, with or without isotopes
        if(deutflag.and..not.o18flag) then
          !read for deut only
          read(15,*,end=998) year(2),day(2),ta_in(2),dp_in(2),   &
                         ua_in(2),sw_in(2),rlwd_in(2),ps_in(2),       &
                         prec_in(2), deutprec_in(2), &
                         runin_in(2), deutrun_in(2)
        else if(o18flag.and..not.deutflag) then
          !read for o18 only
          read(15,*,end=998) year(2),day(2),ta_in(2),dp_in(2),   &
                         ua_in(2),sw_in(2),rlwd_in(2),ps_in(2),       &
                         prec_in(2), o18prec_in(2), &
                         runin_in(2), o18run_in(2)          
        else if(o18flag.and.deutflag) then
          !read for deut and o18      
          read(15,*,end=998) year(2),day(2),ta_in(2),dp_in(2),   &
                         ua_in(2),sw_in(2),rlwd_in(2),ps_in(2),       &
                         prec_in(2), deutprec_in(2), o18prec_in(2), &
                         runin_in(2), deutrun_in(2), o18run_in(2)
        else
          !only wb_flag
           read(15,*,end=998) year(2),day(2),ta_in(2),dp_in(2),   &
                         ua_in(2),sw_in(2),rlwd_in(2),ps_in(2),       &
                         prec_in(2),runin_in(2)
        end if
      else
        !not variable lake depth, no isotopes
        read(15,*,end=998) year(2),day(2),ta_in(2),dp_in(2),   &
                         ua_in(2),sw_in(2),rlwd_in(2),ps_in(2)
        runin_in(2)=0.0
        prec_in(2)=0.0
      end if

!     End Ashling

      if (ispin.le.nspin.and.year(2).ne.year(1)) then
        ispin = ispin+1
        rewind 15
        goto 150
      end if

      call datain (ta_in(2),dp_in(2),ua_in(2),sw_in(2),rlwd_in(2),   &
                   ps_in(2),prec_in(2),runin_in(2),qa_in(2),rh_in(2))

      nsteps = 24*60*60/int(dt) * int((day(2)-day(1)))
!     Begin Ashling
!     Added isotopes_i
      do j=1,nsteps
         xtime = real(j)/nsteps
         call tendency (j,nsteps,ta_in,ta_i,qa_in,qa_i,ua_in,ua_i,   &
                        sw_in,sw_i,rlwd_in,rlwd_i,ps_in,ps_i,        &
                        rh_in,rh_i,prec_in,prec_i,deutprec_in, &
                        deutprec_i,o18prec_in, o18prec_i, &
                        runin_in,runin_i,deutrun_in, deutrun_i,&
                        o18run_in,o18run_i)

         call lake_main(xtime,day(1),ta_i,ua_i,qa_i,ps_i,  &
                        prec_i,sw_i,rlwd_i,runin_i,rh_i,nsteps,&
                        deutprec_i,o18prec_i,deutrun_i,o18run_i)
      enddo
!     End Ashling

      year(1) = year(2)                 ! increment year,day,hour
      day(1) = day(2)

      goto 150
 998  continue

      close(51)
      stop
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    LAKE_MAIN
!    main lake subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !   call lake_main(xtime,day(1),ta_i,ua_i,qa_i,ps_i,  &
      !                  prec_i,sw_i,rlwd_i,runin_i,rh_i,nsteps,&
      !                  deutprec_i,o18prec_i,deutrun_i,o18run_i)
      !Ashling: added isotopes_i (call above)
      subroutine lake_main (xtime, julian, ta_i, ua_i, &
                            qa_i, ps, prec, sw, rlwd, runin,rh,nsteps,&
                            deutprec_i,o18prec_i,deutrun_i,o18run_i)

      implicit none
      include 'Tanganyika.inc'
      real lnetw,lneti,luw,lui,lu,julian,ps,prec,sw,rlwd,      &
           rh,runin,tempice,hice,hsnow,fracice,t,tcutoff,u2w,t2w,     &
           q2w,ti,u2i,t2i,q2i,ta_i,qa_i,ua_i,ta,ua,qa,taC,albs,albi,  &
           albw,swi,sww,tin,tcutk,hicedum,delq,evapw,qhw,qew,evapi,   &
           qhi,qei,delqs,tkw,tki,evap,qe,qh,qbot,qw,qnetice,evaps,    &
           dnsty,de,fracadd,surf_a,fracprv,o18prec,o18run,     &
           deutprec,deutrun,runout,hsprv,snowmelt,o18snow,deutsnow,   &
           tdepth1,tdepth2,xtime,rain,snowadd, &
           deutrun_i, deutprec_i, o18prec_i, o18run_i !Ashling
      integer mixmax,i_shuf,k,iwater,mixdep,j,dep_inc,islice,&
              isave_d,nsteps
      logical snow_flag,melt_flag
      dimension t(max_dep,2), de(max_dep), dnsty(max_dep), ti(max_dep,2)

! =================================================================
!        1. initialize and read in info from previous dt
! =================================================================

      i_shuf = 1   ! grab info from common block
      call shuffle (xtime, julian, i_shuf, depth, d_frac, tempice, &
                    hice, hsnow, fracice, t, snow_flag, o18snow,  &
                    deutsnow, surf_a, ti, mixmax, evap, runout,   &
                    melt_flag, nsteps)

      call salt_init (ps, Tcutoff, trace(1,n_trace)) ! freezing point

      call zero (t2w, t2i, q2w, q2i, u2w, u2i, evapw, evapi, qhw, &
                 qhi, qew, qei, lnetw, lneti, luw, lui, sww, swi)

      if (bndry_flag) then
         if (fracice.lt.1.0)                                           &
            call bndry_flux (ta_i, qa_i, ps, ua_i, t(1,1), Tcutoff,    &
                             u2w, t2w, q2w, hice)
         if (fracice.gt.0.0)                                           &
            call bndry_flux (ta_i, qa_i, ps, ua_i, ti(1,1), Tcutoff,   &
                             u2i, t2i, q2i, hice)
         ta = t2w*(1-fracice) + t2i*(fracice)
         qa = q2w*(1-fracice) + q2i*(fracice)
         ua = u2w*(1-fracice) + u2i*(fracice)
      else
         ta = ta_i
         qa = qa_i
         ua = ua_i
      end if

      taC= ta - 273.15     ! convert weighted air temp to C

      do k=1,depth
        t(k,2) = t(k,1)
        ti(k,1)= t(k,1)
        ti(k,2)= ti(k,1)
      enddo

! ======================================================================
!     2. Calculate added precip and isotopic values- and if snow
! ======================================================================
      hsprv = hsnow
      fracprv = fracice
      rain = prec   ! holds rain (as opposed to snow) amount
      snowmelt = 0.0

! Begin Ashling (no edit yet)
      if (prec.gt.0.0 .and.o18flag) then
         o18prec = o18prec_i!co18prec(taC)
         o18run = o18run_i!co18run(taC)
      end if

      if (prec.gt.0.0 .and.deutflag) then
         deutprec = deutprec_i!cdeutprec(taC)
         deutrun = deutrun_i!cdeutrun(taC)
      end if
! End Ashling

      if (taC.le.snowcut.and.hice.gt.0.0 .and.prec.gt.0) then
        hsnow = hsnow+prec*(rhowat/rhosnow)
        rain = 0.
        snow_flag = .true.
        if (o18flag)  &
           o18snow = (o18snow*hsprv*(rhosnow/rhowat)+prec*o18prec) &
                      /(hsnow*(rhosnow/rhowat))
        if (deutflag) &
           deutsnow = (deutsnow*hsprv*(rhosnow/rhowat)+prec*     &
                       deutprec)/(hsnow*(rhosnow/rhowat))
      endif

!=======================================================================
!      3. Calculate incoming sw radiation over water and ice
!=======================================================================

       call lake_albedo(taC, tcutoff, julian, albs, albi, albw, &
                        melt_flag)
       if (hsnow.gt.snocrit) then
        swi=sw*(1.-albs)    ! swi because over ice fraction
       else if (hsnow.gt.0.0 .and.hsnow.le.snocrit) then
        albi=(albi+albs)/2. ! if thin snow, avg albedos
        swi=sw*(1.-albi)
       else if (hice.gt.0.0 .and.hsnow.le.0.0) then
        swi=sw*(1.-albi)
       endif
       sww=sw*(1.-albw)

!======================================================================
!      4. calculate evaporation, adjust fluxes over ice, adjust for salt
!======================================================================
! ==== 4.1 calculate fluxes for open water ============================

      tin=t(1,1)+273.15     ! surface water T in kelvin
      Tcutk=Tcutoff+273.15  ! freezing point in kelvin
      hicedum=0.0           ! send ice=0.0 to latsens for open water calc
      call latsens (tin, Tcutk, hicedum, ta, qa, ua, ps, delq,  &
                    evapw, Qhw)
      Qew = -evapw*Le

! ==== 4.2 adjust fluxes for ice cover ================================

      if (hice.gt.0.0) then  ! if ice present
        tempice=tempice+273.15
        evapi=evapw
        call adjust_flux (swi, tempice, Tcutk, hice, hsnow, ta, &
                          qa, ua, ps, delq, evapi, Qhi, rlwd)
        Qei = -evapi*Lei ! latent over ice, should be Lei
        tempice=tempice-273.15
      endif ! end if there is ice present

! ==== 4.3 adjust fluxes for salinity  ============================

      if (trace(1,n_trace).gt.0.0) then  ! adjust for salinity
        call salt_evap (trace(1,n_trace), evapw, qa, delq,  &
                        ps, Tcutoff, hice, t(1,1), delqs )
        Qew = -evapw*Le
      endif

!============================================================
!     5. Calculate long wave fluxes over ice and water
!============================================================

      tkw=273.15+t(1,1)
      tki=273.15+tempice
      luw= -0.97*delta*tkw**4. ! long wave up from water surface
      lui= -0.97*delta*tki**4. ! long wave up from ice  surface
      lnetw=rlwd+luw ! net long wave over water
      lneti=rlwd+lui ! net long wave over ice

!=============================================================
!      7. Calculate change in ice thickness and fraction
!         within fraction that already has ice
!=============================================================

      if (fracice.gt.0.0 .or. hsnow .gt. 0.0)  &
        call lake_ice (rlwd, tempice, Qhi, Qei, Tcutoff, swi, hice, &
                       hsnow, snowmelt, ti(1,1), qbot, qw, evapi, &
                       qnetice, fracice, evaps)
!     if (fracice.eq.0.0.and.hsnow.le.0.0) &
      qnetice=0. ! set to zero, used only for iceform now
      fracice=amax1(0.,fracice)
      if (fracice .eq. 0.) hice = 0.0

!==============================================================
!     8.  Adjust temps of water column in open water fraction
!==============================================================

      if (fracprv.lt.1.0) then ! at least some open water
        iwater=1 ! signal that open water calculation
        call eddy (iwater, ua, t, de, depth, trace)
        call temp_profile (iwater, qbot, qw, t, sww, lnetw, Qew, &
                           Qhw, de, depth, trace)
        call tracer_profile (de, depth, iwater)
        mixdep  = 1
        call tracer_mixer(t, dnsty, depth, trace, mixdep, iwater)
        if (mixdep.gt.mixmax) mixmax=mixdep
      endif  ! if there is open water present

!==============================================================
!     9.  Adjust temps of water column in ice fraction
!==============================================================

      if (fracprv.gt.0.0) then ! if  there is ice present
        iwater=0 ! signal that not an open water calculation
        call eddy (iwater, ua, ti, de, depth, trace_i)
        call temp_profile (iwater, qbot, qw, ti, swi, lneti, Qei,  &
                           Qhi, de, depth, trace_i)
        call tracer_profile (de, depth, iwater)
        mixdep = 1
        call tracer_mixer (ti, dnsty, depth, trace_i, mixdep, iwater)
        if (mixdep.gt.mixmax) mixmax = mixdep
      endif  ! if there is ice fraction

!==============================================================
!     10. Calculate ice formation in open water fraction
!==============================================================

      if (fracprv.lt.1.0 .and.t(1,1).lt.Tcutoff) then
        if (iceflag) then
          call ice_form (ps, qnetice, t, depth, Tcutoff, fracprv,  &
                         trace, fracadd, fracice, hice)
         fracice = fracice+fracadd  ! add to frac from lakeice
         hsnow = hsnow*fracprv/fracice  ! conserve snow
        endif
        if (fracadd.eq.-999.) fracice=1.0
      endif  ! if open water for ice to grow in

!  if too much snow for ice buoyancy, switch snow to ice
      snowadd = hice*(rhoice/rhowat-1.)+hsnow*(rhosnow/rhowat)
      if (snowadd.lt.0.) snowadd=0.0
      hice = hice+snowadd
      hsnow = hsnow-snowadd*(rhowat/rhosnow)

!==============================================================
!     11.  Average ice and water columns
!==============================================================

      call column_avg (depth, t, ti, trace, trace_i, fracprv)
      call tracer_avg (depth, fracprv)

!==============================================================
!      12.1 DO WATER AND SALT BALANCE
!==============================================================

      if (wb_flag) then

      tdepth1 = t(depth,1)
      tdepth2 = t(depth,2)
      isave_d = depth

      call water_balance (depth, d_frac, rain, evapw*(1.-fracprv),  &
                          evapi*fracprv, snowmelt, runin, runout,   &
                          surf_a, dep_inc, isave_d)

      if (dep_inc.gt.0) then  ! check about adding new lake layers
         do islice=1,dep_inc  ! add new lake layers
           t(isave_d+islice,1)=tdepth1 ! set temperature of new layers
           t(isave_d+islice,2)=tdepth2
         enddo
      endif

      endif

!==============================================================
!   12.2  DO d18O AND dH BALANCE FOR OPEN WATER ONLY
!==============================================================

      if (o18flag) then
        call O18 (runin, runout, surf_a, rain, evapw*   &
                 (1.-fracprv), rh, t, o18prec, o18run,  &
                 o18snow, snowmelt*rhosnow/rhowat)
      endif
      if (deutflag) then
        call D2H (runin, runout, surf_a, rain, evapw*   &
                 (1.-fracprv), rh, t, deutprec, deutrun,&
                 deutsnow, snowmelt*rhosnow/rhowat)
      endif
      if (hsprv.gt.0.0 .and.hsnow.eq.0) then  ! all snow gone this dt
         snow_flag = .false.  ! no snow present anymore
         o18snow = 0.0
         deutsnow = 0.0
      endif

!==============================================================
!      13.  Place updated info back into common blocks
!==============================================================

      i_shuf = 2   ! put/get info in the common block

      evap = evapw*(1.-fracprv)+evapi*(fracprv)
      if (snowmelt.lt.0) melt_flag=.True.
      if (snowmelt.ge.0) melt_flag=.False.
      call shuffle (xtime,julian,i_shuf, depth, d_frac, tempice,  &
                    hice, hsnow, fracice, t, snow_flag, o18snow,  &
                    deutsnow, surf_a, ti, mixmax, evap, runout,   &
                    melt_flag, nsteps)

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                       ADJUST_FLUX
!   solve for the temperature of the ice (tp) necessary
!   to bring into balance the meteorological heat fluxes and
!   the heat flux in the upper component of the ice/snow.
!   also adjusts sensible heat flux to reflect this new ice
!   temperature(hsen). Evap is changed inasmuch as the surface
!   vapor pressure on which it depends is a function of the new
!   surface temperature.  The subroutine iterates until the solved
!   surface T derived from the interplay of the fluxes is equal
!   to the T fed in to the iteration loop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine adjust_flux (sw, tp, tcutoff, hice, hsnow, ta, qa, &
                              ua, psurf, delq, evap, hsen, rlwd)
      implicit none
      include 'Tanganyika.inc'
      real sw,tp,tcutoff,hice,hsnow,ta,qa,ua,psurf,delq,evap,  &
           hsen,rlwd,t4,x,a,b,qsen,qlat,qmet,condbar,val,val2, &
           q0,t0,tcutc,tposs,t,tlat
      integer ntimes

      t4(x)=(x+273.15)**4.
      tcutc=tcutoff-273.15
      ntimes=1
      a=-0.1
      b=10.0
 99   continue
      tposs=b
199   t=tposs
      tlat=t+273.15   ! switch to kelvin
      call latsens (tlat,tcutoff,hice,ta,qa,ua,psurf,delq,evap,qsen)
      hsen=qsen
      qlat=-evap*Lei
      qmet=rlwd-emis*delta*t4(t)+qsen+qlat
      call ice_rad (sw, hice, hsnow, condbar, val, val2)
      q0=-qmet
      t0=condbar*(sw-q0)+tcutc-val ! eq 7 from p + h
      if (t0.ge.t.and.ntimes.eq.1) then
        b=t-a
        a=-.001
        ntimes=2
        goto 99
      else if (t0.ge.t.and.ntimes.eq.2) then
        b=t-a
        a=-.00001
        ntimes=3
        goto 99
      else if (t0.ge.t.and.ntimes.eq.3) then
        tp=t+273.15         ! temperature to pass back out
        goto 299
      endif
      tposs=tposs+a
      if (tposs.lt.-100.) goto 299
      goto 199
 299  continue

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                         BNDRY_FLUX
!     used to interpolate lowest atmospheric model level to 2-meters
!     based on GENESIS/LSX codes
!     the flag in the include file determines when this is actually used
!     the computations are based on the sigma level of the gcm as
!     specified in include file and should only be called when using
!     atmospheric model output
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine bndry_flux(ta,qa,psurf,ua,tlake,Tcutoff,u2,t2,q2,hice)

      implicit none
      include 'Tanganyika.inc'
      real ta,qa,psurf,ua,tlake,Tcutoff,u2,t2,q2,pa,tlakeK,cappa,tfac, &
           cdmaxa,rhoa,cp,aloga,alogg,a,b,eo,qo,richu,stramu,strahu,  &
           cdh,cdm,fsena,fvapa,tau,za,zb,zt,ugrad,tgrad,qgrad,cdmaxb, &
           z,hice
      integer i
      parameter (cdmaxb=1.e20)

      tlakeK = tlake + 273.15

      cappa = rair/cpair
      tfac = 1./(sigma**cappa)   ! correction factor for temp at za
      cdmaxa = 100./(2.*dt)

      pa = sigma*psurf                    ! pressure at za
      rhoa = pa / (rair*ta*(1.+(rvap/rair-1.)*qa))  ! air density
      cp = cpair*(1.+(cvap/cpair-1.)*qa)    ! specific heat
      za = (psurf-pa)/(rhoa*grav)         ! height of lowest sigma level

      aloga = log(za)
      alogg = log(zo)

      if (hice.le.0.0 .and.tlake.gt.Tcutoff) then  ! for eo and qo calcs
         a=c72
         b=c73
      else          ! if ice, calc eo over ice
         a=c70
         b=c71
      endif

      eo=ca*exp(a*(tlakeK-cb)/(tlakeK-b)) ! surface vapor pressure
      qo=0.622*eo/(psurf-(1.-.622)*eo)  ! surface specific humidity

      call fstrat (tlakeK,ta,tfac,za,alogg,aloga,ua,richu, &
                   stramu,strahu)  ! determine stability

      cdh = ua*(kv/(aloga-alogg))**2*strahu ! calc air transfer coeffs
      cdm = ua*(kv/(aloga-alogg))**2*stramu

      cdh = min (cdmaxa,cdh/(1.+cdh/cdmaxb))
      cdm = min (cdmaxa,cdm/(1.+cdm/cdmaxb))

      fsena = rhoa*cdh*cp*(ta*tfac-tlakeK) ! calc fluxes
      fvapa = rhoa*cdh*(qa-qo)
      tau = rhoa*cdm*ua

      zb=zo
      u2=0.0    ! initialize surface values
      t2=tlakeK
      q2=qo

      do i=1,2
         zt=float(i)
         if (i.eq.2) zb=1.0
         z=(zt-zb)/(log(zt/zb))
         ugrad=(sqrt(tau/rhoa))/(kv*z*sqrt(stramu))
         tgrad=-fsena*(log(zt/zb))**2/(rhoa*ugrad*kv**2*strahu*cp)
         qgrad=-fvapa*(log(zt/zb))**2/(rhoa*ugrad*kv**2*strahu)
         u2=u2+ugrad
         t2=t2-tgrad
         q2=q2-qgrad
      enddo

222   continue

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 COLUMN_AVG
!     average water and ice columns for fractional cover
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine column_avg (depth,t,ti,salty,saltyi,fracprv)

      implicit none
      include 'Tanganyika.inc'
      integer j
      real t(max_dep,2),ti(max_dep,2),dnsty(max_dep)
      real salty (max_dep,n_trace), saltyi(max_dep,n_trace)
      real salt_mx (max_dep),temp,fracprv,dnstyw,dnstyi,cpw,cpi,z

      do j=1,depth
        call density(t(j,1),salty(j,n_trace),dnstyw)
        call density(ti(j,1),salty(j,n_trace),dnstyi)
        call specheat(t(j,1),salty(j,n_trace),cpw)
        call specheat(ti(j,1),saltyi(j,n_trace),cpi)
        dnstyw=dnstyw+1000.
        dnstyi=dnstyi+1000.
        z=dz
        if (j.eq.1) z=surf
        temp=((1.-fracprv)*t(j,1)*z*dnstyw*cpw+   &
          fracprv*ti(j,1)*z*dnstyi*cpi)/          &
          ((z*dnstyw*cpw+z*dnstyi*cpi)*0.5)

        salt_mx(j) = (1.-fracprv)*salty(j,n_trace) +  &
                 fracprv*saltyi(j,n_trace)
        t(j,1)=temp
        call density (t(j,1),salt_mx(j),dnsty(j))
      enddo

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                        D2H
!     computes deuterium fractionation from surface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       subroutine D2H(runin, runout, surf_a, rain,   &
                      evap, rh, t, deutprec,deutrun, &
                      deutsnow,hs)
       implicit none
       include 'Tanganyika.inc'
       real runin,runout,surf_a,rain,evap,rh,deutrun,runout_len,  &
            runindeut,runoutdeut,deutprec,tk,alphadeut,rldeut,    &
            radeut,redeut,delvdeut,evapdeut,runin_len,precdeut,   &
            hs,snowdeut,deutsnow
       real t(max_dep,2)

       tk=t(1,1)+273.15
       alphadeut = exp(24844./(tk**2.)- 76.248/tk + 0.05261) ! activity coef
       rldeut = (trace(1,2)*1.e-03) + 1.
       radeut = (deutair*1.e-03) + 1.
       redeut=(rldeut/alphadeut-rh*f*radeut) / (((1.-rh)/alphak) + &
              (rh*(1.0-f)))
!      redeut is the fractionation factor for isotopes in boundary layer air
       delvdeut = (redeut-1.)*1.e03
       evapdeut = delvdeut * evap
       runoutdeut  = runout * trace(1,2)
       snowdeut= hs*rhosnow/rhowat*deutsnow
       runin_len = (runin *(b_area-surf_a) / surf_a)
       runindeut = runin_len * deutrun
       precdeut = rain * deutprec
       trace(1,2)=(trace (1,2)*dexch+   &
          runindeut-runoutdeut+precdeut-evapdeut/1.e3*dt+snowdeut)/ &
          (dexch + runin_len-runout + rain - evap/1.e3*dt+hs*rhosnow/&
           rhowat)

       return
       end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                   DATA_IN
!  convert input data into form needed for model
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine datain (ta_in,dp_in,ua_in,sw_in,rlwd_in,ps_in, &
                         prec_in,runin_in,qa_in,rh_in)
      implicit none
      include 'Tanganyika.inc'
      real ta_in,qa_in,dp_in,ua_in,ps_in,prec_in,runin_in,sw_in,   &
           rlwd_in,es,ea,rh_in

      if(sw_in.lt.0.0)sw_in=0.0
!      ps_in = ps_in*100.                      ! convert from mb to Pa   
      prec_in = prec_in/1000.                 ! convert from mm to m
      runin_in =runin_in/1000.               ! convert from mm to m

!      es = ca*exp(c72*(ta_in+273.16-cb)/(ta_in+273.16-c73))   ! sat vap
      es = ca*exp(c72*(ta_in-cb)/(ta_in-c73))   ! sat vap
!     ea = ca*exp(c72*(dp_in+273.16-cb)/(dp_in+273.16-c73))   ! act vap pres
!     qa_in=ea*0.622/(ps_in-ea*0.378)
!      rh_in = ea/es
      rh_in = dp_in
      rh_in = rh_in/100.
      ea = rh_in*es
!      qa_in = dp_in
!      rh_in = -999.

      qa_in = ea*0.622/(ps_in-ea*0.378)
!      ta_in = ta_in+273.15
      ua_in = ua_in*log(z_screen/0.001)/log(10.0/0.001)  ! use for 10-meter wind
!      write (*,*) dp_in, es, ea, rh_in, qa_in
      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  DENSITY
! calculate density as a function of temperature and salinity
! temp is deg C. salinity (s) is ppt, and rhostp is in kg /m3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine density (ts,ss,rhostps)
      implicit none
      real*8 rhow, rhost, t, s, rhostp
      real ts,ss,rhostps
      t=dble(ts)
      s=dble(ss)

      rhow=999.842594D0+6.793952D-2*t-9.095290D-3*t**2   &
          +1.001685D-4*t**3-1.120083D-6*t**4+6.536332D-9*t**5

      rhost=rhow+s*(0.824493D0-4.0899D-3*t+7.6438D-5*t**2  &
         -8.2467D-7*t**3+5.3875D-9*t**4)+                  &
           s**(3.D0/2.D0)*(-5.72466D-3+1.0227D-4*t         &
          -1.6546D-6*t**2)+4.8314D-4*s**2

      rhostp=rhost  ! no pressure version
      rhostps=rhostp-1.D3

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                     EDDY
!   computes eddy diffusion profile
!   calculation from Henderson-Sellers 1985 Applied
!   Mathematical Modelling
!   cm corrected code for rad calculation in order to match
!   Henderson-Sellers publication 2/4/2016
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine eddy (iwater,u2,t,de,depth,salty)

      implicit none
      include 'Tanganyika.inc'
      real u2,u,ks,N2,ws,Po,radmax,z,rad,Ri,dpdz,zhalf
      real de(depth), dnsty(depth), t(depth,2)
      real salty(max_dep,n_trace)
      integer k,iwater

      do k=1,depth
         call density (t(k,1),salty(k,n_trace),dnsty(k))
      enddo

      if (iwater.ne.1) then ! if ice fraction, no de just dm
        do k=1,depth
          de(k)=dm
        enddo
        return !  no further calculations needed
      endif

      u=amax1(u2,0.5) ! avoid NAN in ks
      ks=6.6*sqrt(abs(sin(xlat*raddeg)))*u**(-1.84)
      ws=0.0012*u2
      Po=1.0

      radmax=4.e4 ! limits Ri to 10
      do k= 1,depth-1
        if (k .eq. 1) then  ! for surface layer different depth
          zhalf=(surf+dz)*0.5
        else
          zhalf=dz
        endif
        dpdz=(dnsty(k+1)-dnsty(k))/zhalf
        N2=(dpdz/(1.e3+dnsty(k)))*grav
        z=surf+float(k-1)*dz
        if ((ks*z)/ws.gt.40.) then
          rad = radmax  ! avoid NAN
        else
          rad=1.+40.*N2*(kv*z)**2./(ws**2.*exp(-2.*ks*z))
          if (rad.gt.radmax) rad=radmax
        endif
        rad=amax1(rad,1.0) ! so that Ri lower lim is 0.0
        Ri=(-1.0+sqrt(rad))/20.0
        if ((ks*z)/ws.gt.40.) then
           de(k)=dm
        else
           de(k)=dm+kv*ws*z*Po*exp(-ks*z)/(1.0+37.0*Ri**2)
        end if
      enddo
      de(depth)=dm
      de(depth-1)=de(depth)  ! necessary for cn solution to work

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 FILE_OPEN
!  open input and output files
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine file_open
      !Begin Ashling
      implicit none
      include 'Tanganyika.inc'

!      input files
      open(unit=15,file=datafile, status='old')
      !End Ashling

!      output files
      open(unit=50,file='ERA-HIST-Tlake_prof.dat',status='unknown')
      open(unit=51,file='ERA-HIST-Tlake_surf.dat',status='unknown')

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                    FSTRAT
!    computes mixing-length stratification correction factors
!    for momentum and heat/vapor, based on genesis lsx code
!    used by BNDRY_FLUX for interpolation from atmo model
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine fstrat(tb,tt,ttfac,zt,alb,alt,u,rich,stram,strah)

!     tb=bottom temperature
!     tt=top temperature
!     ttfac=potential temperature correction to tt
!     zb=height of bottom
!     zt=height of top
!     alb=log(bottom roughness length)
!     alt=log(z at top)
!     u=wind speed at top
!     rich=richardson number
!     stram=stratification factor for momentum
!     strah=stratification factor for heat/vap
      implicit none
      include 'Tanganyika.inc'

      real tb,tt,ttfac,zt,alb,alt,u,rich,stram,strah,zb,x,c,sqri

      zb = zo
      rich=grav*max(zt-zb,0.)*(tt*ttfac-tb)/(tt*ttfac*u**2)
      rich=min(rich,1.0)
      if (rich.le.0) then
         x=max(alt-alb,0.5)
         c=(kv/x)**2*9.4*exp(0.5*x)
         sqri=sqrt(-rich)
         stram=1.-9.4*rich/(1.+7.4*c*sqri)
         strah=(1.-9.4*rich/(1.+5.3*c*sqri))/0.74
      else
         stram=1./(1.+4.7*rich)**2
         strah=stram/0.74
      endif

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 ICE_FORM
! calculates the fractional coverage of new ice that
! has formed over the open water during this timestep
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine ice_form (psurf,qnetice,t,depth,Tcutoff,  &
                          fracprv,salty,fracadd,fracice,hi)
      implicit none
      include 'Tanganyika.inc'
      real t(max_dep,2), salty(max_dep,n_trace)
      real qnetice,Tcutoff,fracprv,fracadd,fracice,hi,sum,extra,  &
           cp,di,xfrac,dnsty,psurf
      integer j

      call salt_init (psurf,Tcutoff,salty(1,n_trace))
      sum=0.
      do j=1,depth
        if (t(j,1).lt.Tcutoff) then
           call density (t(j,1),salty(j,n_trace),dnsty)
           call specheat (t(j,1),salty(j,n_trace),cp)
           extra=(Tcutoff-t(j,1))*dz*(dnsty+1.e3)*cp
           if (j.eq.1) extra=(Tcutoff-t(j,1))*surf*(dnsty+1.e3)*cp
           t(j,1)=Tcutoff
           sum=sum+extra
        endif
      enddo
      qnetice=(sum/dt)*(1.0-fracprv) ! heat flux absorbed into ice
      if (fracprv.le.0.0) hi=fracmin
      di=sum/(fusion*rhoice) ! thickness of new ice
      fracadd=(di/fracmin)*(1.0-fracprv) ! convert to fracadd
      if ((fracadd+fracice).gt.1.0) then  ! too much added
        xfrac=(fracice+fracadd)-1.0
        di=xfrac*fracmin/1.0
        hi=hi+di
        fracadd=-999.
      endif

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                      ICE_RAD
!  calculates terms needed for surface energy balance in
!  presence of lake ice and snow
!  equations based on Patterson and Hamblin (1988) Limnology
!  and Oceanography v. 33 p. 323
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine ice_rad (sw, hi, hs, condbar, val, val2)

      implicit none
      include 'Tanganyika.inc'
      real  a,b,c,d,sw,hi,hs,val,val2,condbar

      condbar=(hs*condi+hi*conds)/(condi*conds)
      a=(1.-exp(-lamssw*hs))/(conds*lamssw)
      b=exp(-lamssw*hs)*(1-exp(-lamisw*hi))/(condi*lamisw)
      c=(1.-exp(-lamslw*hs))/(conds*lamslw)
      d=exp(-lamslw*hs)*(1-exp(-lamilw*hi))/(condi*lamilw)
      val=sw*afrac1*(a+b)+sw*afrac2*(c+d)
      val2=-afrac1*sw*(1-exp(-(lamssw*hs+lamisw*hi))) &
            -afrac2*sw*(1-exp(-(lamslw*hs+lamilw*hi)))

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                INIT_LAKE
!    initialize lake
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine init_lake

      implicit none
      include 'Tanganyika.inc'
      integer k

!==============================================================
!     1. Initialize depth and salinity
!==============================================================

      depth_a = depth_begin
      salty_a = salty_begin

!==============================================================
!     2. Initialize lake temp and tracer profiles
!===============================================================

! Begin Ashling
      do k=1,depth_a
          temp_a(k)= tempinit 
          trace_a(k,1)= o18init 
          trace_a(k,2)= deutinit 
          trace_a(k,3)= salty_a  ! salinity
      enddo
 ! End Ashling

      if (s_flag) then ! salinity varies
         num_tra = n_trace  ! do tracer calculations for all
      else
         num_tra = n_trace - 1 ! don't do salinity
      endif

!==============================================================
!     3. Zero others
!==============================================================

      d_fraca = 0.   ! incremental lake level change
      hice_a = 0.    ! height of lake ice
      hsnow_a = 0.    ! height of snow on lake ice
      fraci_a = 0.    ! fraction of lake covered by ice
      tempi_a = 0.    ! temperature of lake ice
      snow_flag_a = .false.  ! snow present
      melt_flag_a = .false.
      o18snow_a = 0.0    ! accumulated d18O in snow on lake ice
      deutsnow_a = 0.0   ! accumulated dD in snow on lake ice
      surfarea_a = 0.0  ! surface area of lake

      mix_ave = 0.0  ! holds mixing depth for time ave
      tsurf_ave = 0.0  ! holds lake surface T for time ave
      fice_ave = 0.0   ! holds ice fraction for time ave
      evap_ave = 0.0   ! holds lake evap for time ave
      hice_ave = 0.0   ! holds ice height for time ave

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  LAKE_ALBEDO
!  From BATS, calculates fragmented albedos (direct and diffuse) in
!  wavelength regions split at 0.7um.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine lake_albedo (tair, tcutoff, julian, albs, albi, albw, &
                              melt_flag)

      implicit none
      include 'Tanganyika.inc'
      real tair,tcutoff,albs,albi,albw,tdiffs,albgl,   &
           albgs,tdiff,julian
      logical melt_flag

      tdiffs=tair-tcutoff  ! degrees above freezing in celcius
      tdiff=amax1(tdiffs,0.) ! make sure tdiffs above zero
      tdiffs=amin1(tdiff,20.) ! limit diff to be < 20 degrees
      albgl=sical1-1.1e-2*tdiffs  ! long wave = near-infrared
      albgs=sical0-2.45e-2*tdiffs ! short wave = visible
      albi=fsol1*albgs+fsol2*albgl ! wt. long.v.short by fsol

!     snow albedo really depends on snow-age, zenith angle,
!     and thickness of snow, age reduces visible rad snow albedo

      if (melt_flag) then
         albs=alb_slush
      else
         albs=alb_snow  ! albedo of snow
      end if

      albw = 0.08+0.02*sin(2.*pi*julian/365.+pi/2.)

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                    LAKE_DRAG
! calculates drag coefficient using BATS method
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine lake_drag (tlakek, t1k, u1, cdrx)

      implicit none
      include 'Tanganyika.inc'
      real tlakek,t1k,u1,cdrx,cdrmin,cdr,ribn,ribd,rib,zs

!     cdrn=(kv/alog(z2/zo))**2.  ! neutral drag coefficient
      zs = z_screen
      if (bndry_flag) zs = 2.0
      ribn=zs*grav*(1.-tlakek/t1k)  ! neutral bulk Richardson number

      if ((tlakek/t1k).le.1.0) then  ! unstable conditions
        ribd=u1**2.+0.1**2.  ! other bulk Richardson number
      else  ! stable conditions
        ribd=u1**2.+1.0**2.
      endif

      rib=ribn/ribd   !  bulk Richardson number
      if (rib.lt.0.) then ! unstable conditions
        cdr=cdrn*(1.0+24.5*(-cdrn*rib)**0.5)
      else   ! stable conditions
        cdr=cdrn/(1.0+11.5*rib)
      endif
      cdrmin=amax1((.25*cdrn),6.e-4)

      if (cdr.lt.cdrmin) cdr=cdrmin
      cdrx=cdr

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                LAKE_ICE
!  adds to or subtracts from lake ice and snow thickness
!  based on evaporation /condensation and melting/freezing
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine lake_ice (radlwd,tempice,qsen,qlat,      &
                           tcutoff,sw,hi,hs,ds,twater,qbot,  &
                           qw,evapi,qnetice,fracice,evaps)
      implicit none
      include 'Tanganyika.inc'
      real radlwd,qsen,qlat,tcutoff,sw,hi,hs,twater,qbot,      &
           qw,ds,evapi,qnetice,fracice,evaps,tmelt,t4,x, &
           q0t0,condqw,evapl,tprev,qmet,val,val2,q0,tempice,   &
           qmelts,condbar,qf,qmeltb,qmeltsx,dibot,hiprv,         &
           extradi,df,xfrac,di,diextra,extraf,disurf
      parameter (tmelt = 0.0)  ! melting temp for snow

      t4(x)=(x+273.15)**4.
      q0t0(x)=sw+(1./condbar)*(tcutoff-(x)-val)
!     conductivity for water to ice flux
!=============================================================
!     1.  Calculate surface fluxes, update tempice
!=============================================================

       condqw = rhowat * cpw_ice * surf / (qwtau*2.)
       evapl=evapi*dt/1000.  ! convert from mm/sec to m (over this dt)
       tprev=tempice          ! keep track of incoming t0
       qmet=radlwd-emis*delta*t4(tprev)+qsen+qlat
!     1a. Calculate surface fluxes, update tempice
!         qmet is the same as H(To) in Patterson and Hamblin, pg. 327, eqn. 6
       call ice_rad (sw, hi, hs, condbar, val, val2 )
       q0=-qmet  ! q0 is flux of heat out of i/s from i/s surf
       tempice=condbar*(sw-q0)+tcutoff-val ! eq 7 from p + h
       qbot=sw+val2   ! this is sw rad leaving base of ice

!=============================================================
!     2.  Adjust tempice if greater than melt temp
!=============================================================

       if (tempice.gt.tmelt) then
         q0=q0t0(tmelt)     ! set q0 to q0 for t0 = tmelt
         tempice=tmelt      ! set surf temp to melting temp
         qmelts=q0+qmet     ! calc extra heat flux for melting
       else
         qmelts=0.0
       endif

!==============================================================
!     3.  Calculate fluxes at the base of the ice
!==============================================================

       qf=q0+val2  ! flux of heat into ice from base of ice
       qw=-condqw*(tcutoff-twater)
       qmeltb=qf-qw  !  note opposite sign from surface qmelt
       qnetice=qmeltb-qmelts  ! flux for freeze/melt

!==============================================================
!     4.1  Adjust snow depth for evap/condensation
!==============================================================

      qmeltsx=0.0           ! initialize here, incase no snow
      if (hs.gt.0.0) then
        if (evapl*(rhowat/rhosnow).le.hs) then
          hs=hs-evapl*(rhowat/rhosnow)   ! not water equivalent
          evapl=0.0                   ! all evapl used in removing snow
          evaps=evapi                 ! all ice evap to snow
        else     ! evapl is greater than snow thickness
          evaps=hs*(rhosnow/rhowat)*1.e3/dt  ! make a rate
          evapl=evapl-hs*(rhosnow/rhowat)
          hs=0.0
        endif
      else
        evaps = 0.0
      endif

!===============================================================
!    4.2  Adjust snow depth for melting/freezing
!===============================================================

      if (hs.gt.0.0) then
        ds=(-qmelts/(rhosnow*fusion))*dt  ! ds < 0, melting
        if (-ds.gt.hs) then  ! then have to melt ice too
          qmeltsx=qmelts-(hs*rhosnow*fusion/dt)  ! enery remaining for ice
          hs=0.0             ! set snow to zero
        else
          hs=hs+ds
        endif                ! -ds > hs
      endif                 ! if hs gt 0.0

!===============================================================
!    4.3   Calculate ice thickness change for m/f and e/c
!===============================================================

      if (hs.le.0.0) then   ! if there is ice at the surface
        disurf=((-qmelts/(rhoice*fusion))+     &
               (-qmeltsx/(rhoice*fusion)))*dt+ & ! add extra from snow
               (-evapl*(rhowat/rhoice))      ! add remaining evap length
      else
        disurf=0.0
      endif
      dibot=(qmeltb/(rhoice*fusion))*dt  ! no minus here

!==================================================================
!     4.4  Adjust ice thickness and fraction
!==================================================================

        if (fracice.ge.1.0) then
          hiprv=hi
          hi=hi+disurf+dibot ! if full ice, change thickness
          if (hi.lt.fracmin) then
            extradi=fracmin-hi
            df=(extradi/fracmin)*fracice  ! fracice == 1. here
            fracice=1.0-df  ! change fraction cause down to min thickenss
            hi=fracmin  ! set to min thickness
          endif
        else    ! fractional ice cover
          df=fracice*(disurf+dibot)/fracmin ! OK2use fracice, cause=fprv
          fracice=fracice+df
          if (fracice.gt.1.0) then ! too much area, make thicker
            xfrac=fracice-1.0
            di=xfrac*fracmin/1.0
            hi=hi+di
            fracice=1.0              !set frac to 100%
          endif
          if (fracice.lt.fraclim.and.df.le.0.0) then
            xfrac=fracice ! all remaining ice is extra
            diextra=xfrac*fracmin/1.0 ! convert this to a thickness
            extraf=-diextra*rhoice*fusion*(1./dt)
            qw=qw-extraf  ! adjust flux from water for heating
            qnetice=qnetice+extraf
            fracice=0.0  ! set frac, thickness and snow to 0.0
            hi=0.0
            ds = ds+hs   ! remaining snow goes into lake
            hs=0.0       ! no more snow left on ice
          endif
       endif

       evapi=evapl

       return
       end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             LATSENS
!     computes latent and sensible heat fluxes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine latsens (tsurf, Tcutoff, hice, t2,    &
                          q2, u2, psurf, delq, evap, qsen)
      implicit none
      include 'Tanganyika.inc'
      real tsurf,tcutoff,hice,t2,q2,u2,psurf,delq,evap,qsen,  &
           elake,qlake,rai,cdrx,a,b,delt,rash,relhum,pv,pd,rhosurf

      call lake_drag (tsurf, t2, u2, cdrx)   ! calculate drag coef.

      if (hice.le.0.0 .and.tsurf.gt.Tcutoff) then
        a=c72
        b=c73
      else          ! if ice, calc eog over ice
        a=c70
        b=c71
      endif

      elake=ca*exp(a*(tsurf-cb)/(tsurf-b)) ! surface vapor pressure
      qlake=0.622*(elake/(psurf-0.378*elake)) ! surface specific humidity
      relhum=100.*q2/qlake
      pv = (relhum*elake)/100.
      pd = psurf-pv
      rhosurf = pd/(rair*t2)+pv/(rvap*t2)

      delq=q2-qlake  ! specific humidity gradient
      rai=cdrx*u2*rhosurf
      evap=-rai*delq

      delt=t2-tsurf ! temperature gradient
      rash=cdrx*u2*rhosurf*cpair
      qsen=rash*delt

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                        O18
!     computes O18 balance in the surface layer
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine O18 (runin,runout,surf_a,rain,evap,rh,t,  &
                      o18prec,o18run,o18snow,hs)
      implicit none
      include 'Tanganyika.inc'
      real runin,runout,surf_a,rain,evap,rh,tk,alphao18,rlo18,    &
           rao18,reo18,delvo18,evapo18,runin_len,runin_o18,       &
           preco18,o18prec,o18run,runout_len,runout_o18,o18snow,  &
           hs,snowo18
      real t(max_dep,2)
      integer k

      tk=t(1,1)+273.15
      alphao18 = exp(1137./(tk**2.)- 0.4156/tk - 0.00207) ! activity coef
      rlo18 = (trace(1,1)*1.e-03) + 1.
      rao18 = (o18air*1.e-03) + 1.
      reo18=(rlo18/alphao18-rh*f*rao18) / (((1.-rh)   &
            /alphak) + (rh*(1.0-f)))
!     reo18 is the fractionation factor for isotopes in boundary layer air
      delvo18 = (reo18-1.)*1.e03
      evapo18 = delvo18 * evap
      runout_o18 = runout * trace(1,1)
      snowo18 = hs*rhosnow/rhowat*o18snow
      runin_len = (runin*(b_area-surf_a) / surf_a)
      runin_o18  = runin_len * o18run
      preco18 = rain * o18prec
      trace(1,1)=(trace (1,1)*dexch+   &
           runin_o18-runout_o18+preco18-evapo18/1.e03*dt+  &
           snowo18)/ (dexch+runin_len-runout+rain-       &
           evap/1.e03*dt+hs*rhosnow/rhowat)

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 SALT_EVAP
!   calculates change in surface vapor pressure (and
!   therefore evaporation) due to salinity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine salt_evap (salty, evap, qa, delq,   &
                           psurf, Tcutoff, hice, tg, delqs)
      implicit none
      real salty,evap,qa,delq,psurf,Tcutoff,hice,tg,    &
           delqs,qsurf,psurfpa,esurf,elower,esurfs,qgs,alpha

      if (hice.eq.0.0 .and.tg.gt.Tcutoff) then
       qsurf=qa-delq
       psurfpa=psurf
       esurf=(psurfpa*qsurf)/(0.622+0.378*qsurf)
       elower=133.3224*(exp((1.186094*alog(salty))-  &
           (5580.475512/(tg+273.15))+13.674717))
       esurfs=esurf-elower
       qgs=0.622*(esurfs/(psurfpa-0.378*esurfs))
       delqs=qa-qgs
       if (delq.eq.0) delq=1.e-15
         alpha=evap/delq
         evap=alpha*delqs
      endif     ! if ice, no need to lower for salty,

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                   SALT_INIT
! calculates freezing point, in Gill 1982, from Millero 1978
! tfsp is in degrees c, pressure is in decibars
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine salt_init (ps,tfsp,s)

      implicit none
      real tfsp,s,ps,psd

      psd = ps/10000.             ! convert from Pa to decibars
      tfsp=-0.0575*s+1.710523e-3*s**(3./2.)-2.154996e-4*s**2.     &
         -7.53e-4*psd

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                       SHUFFLE
! shuffles data to and from common block
!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine shuffle (xtime,day,i_shuf,depth,d_frac,tempice,   &
                          hice,hsnow,fracice,t_shuf,snow_flag, &
                          o18snow,deutsnow,surf_a,ti_shuf,mixmax,   &
                          evap,runout,melt_flag,nsteps)
      implicit none
      include 'Tanganyika.inc'
      real tempice,hice,hsnow,fracice,surf_a,evap,xtime,day, &
           t_shuf (max_dep,2),ti_shuf(max_dep,2),econv,runout,&
           o18snow,deutsnow
      logical snow_flag,melt_flag
      integer k,i_shuf,i_trace,mixmax,nsteps

      if (i_shuf.eq.1) then  ! get info from common block
        depth   = depth_a
        d_frac  = d_fraca
        tempice = tempi_a
        hice    = hice_a
        hsnow   = amax1(0.,hsnow_a)
        snow_flag = snow_flag_a
        melt_flag = melt_flag_a
        o18snow = o18snow_a
        deutsnow = deutsnow_a
        fracice = fraci_a
        surf_a =surfarea_a
        mixmax = mixmax_a
        do k = 1, depth
           t_shuf(k,1) = temp_a (k)
           ti_shuf(k,1) = ti_a (k)
           do i_trace = 1,n_trace
             trace(k,i_trace) = trace_a (k,i_trace) ! for open water frac
             trace_i(k,i_trace) = trace_a (k,i_trace)  ! for ice frac
           enddo
        enddo
      else if (i_shuf.eq.2) then ! place info in common blocks
        depth_a = depth
        d_fraca = d_frac
        tempi_a = tempice
        hice_a  = hice
        hsnow_a = amax1(0.,hsnow)
        snow_flag_a = snow_flag
        melt_flag_a = melt_flag
        fraci_a  = fracice
        o18snow_a = o18snow
        deutsnow_a = deutsnow
        surfarea_a =surf_a
        mixmax_a =mixmax
        do k = 1, depth
           temp_a (k) = t_shuf(k,1)
           ti_a(k) = ti_shuf(k,1)
           do i_trace = 1,n_trace
             trace_a(k,i_trace) = trace(k,i_trace)
           enddo
        enddo

        mix_ave = mix_ave+mixmax
        tsurf_ave = tsurf_ave+t_shuf(1,1)
        fice_ave = fice_ave+fracice
        evap_ave = evap_ave+evap
        hice_ave = hice_ave+hice
        hsnow_ave = hsnow_ave+hsnow
        o18_ave = o18_ave+trace_a(1,1)
        deut_ave = deut_ave+trace_a(1,2)
        runout_sum = runout_sum+runout
        do k=1,depth
           temp_ave(k)=temp_ave(k)+temp_a(k)
        end do

        if (xtime.eq.1.) then    ! print daily averages
           econv = 60.*60.*24      ! convert from mm/s to mm/day

! Begin Ashling          
! Fixed so that oxygen and hydrogen isotopes can be written separetely 
! Moved format inside of if so that it formats correctly
! I removed ice from my output so that would need to be added back in and the format adjusted accordingly
           if (o18flag.and.deutflag) then
!          write(51,361) day,tsurf_ave/nsteps,&
!           				 fice_ave/nsteps,&
!                         evap_ave*econv/nsteps,&
!                         mix_ave/nsteps,&
!                         hice_ave/nsteps,hsnow_ave/nsteps,&
!                         o18_ave/nsteps, &
!                         deut_ave/nsteps,runout_sum,mixmax,&
!                         depth,d_frac
         write(51,361) day,tsurf_ave/nsteps,&
                         evap_ave*econv/nsteps,&
                         mix_ave/nsteps,&
                         o18_ave/nsteps, &
                         deut_ave/nsteps,runout_sum,mixmax,&
                         depth,d_frac
  361      format(F9.2,1x,6(F6.2,1x),2(I5,1x),F6.2)                                     
!          write(50,*) day,((temp_ave(k)/48.),k=1,depth)
           write(50,362) day,((temp_ave(k)/nsteps),k=1,depth,10)
  362      format(58(F11.2,1x))           
           else if (o18flag.and..not.deutflag) then 
         write(51,371) day,tsurf_ave/nsteps,&
                         evap_ave*econv/nsteps,&
                         mix_ave/nsteps,&
                         o18_ave/nsteps, &
                         runout_sum,mixmax,&
                         depth,d_frac
  371      format(F9.2,1x,5(F6.2,1x),2(I5,1x),F6.2) 
           write(50,372) day,((temp_ave(k)/nsteps),k=1,depth,10)    
  372      format(58(F11.2,1x))                                
           else if (.not.o18flag.and.deutflag) then
         write(51,381) day,tsurf_ave/nsteps,&
                         evap_ave*econv/nsteps,&
                         mix_ave/nsteps,&
                         deut_ave/nsteps, &
                         runout_sum,mixmax,&
                         depth,d_frac
  381      format(F9.2,1x,5(F6.2,1x),2(I5,1x),F6.2)  
           write(50,382) day,((temp_ave(k)/nsteps),k=1,depth,10)   
  382      format(58(F11.2,1x))                               
           else
           write(51,391) day,tsurf_ave/nsteps,mix_ave/nsteps,&
                         evap_ave*econv/nsteps,mixmax,depth
           write(50,392) day,((temp_ave(k)/nsteps),k=1,depth,10)
  391      format(4(F11.4,1x),2(i3,1x))
  392      format(58(F11.2,1x)) !Here change the temperature profile format.           
           end if

! End Ashling

           mix_ave=0.0
           tsurf_ave=0.0
           fice_ave = 0.0
           evap_ave = 0.0
           hice_ave = 0.0
           hsnow_ave = 0.0
           o18_ave = 0.0
           deut_ave = 0.0
           mixmax_a = 1
           runout_sum = 0.0
           do k=1,max_dep
             temp_ave(k)=0.0
           end do
        endif
      endif

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                        SOLAR_DEC
!  computes the solar declination angle from the julian date.
!  from bats version 3 by Dickinson et al.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine solar_dec (julian, declin)

      implicit none
      include 'Tanganyika.inc'
      real declin,obecl,sinob,olong,arg,julian

      obecl=oblq*raddeg
      sinob=sin(obecl)
      if (julian .ge. 81.) olong=dpd*(julian-81.)
      if (julian .lt. 81.) olong=dpd*(julian+284.)
      olong=olong*raddeg  ! longitude from vernal equinox
      arg=sinob*sin(olong)
      declin=asin(arg)   ! solar declination angle

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                       SPECHEAT
! calculates specific heat as a function of salinity and temperature
! t is in deg C, s is salinity ppt, and cpts is J / Kg K
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine specheat (t,s,cpts)

      implicit none
      real t,s,cpts,cpt

      cpt=4217.4-3.720283*t+0.1412855*t**2-2.654387e-3*t**3  &
          +2.093236e-5*t**4
      cpts=cpt+s*(-7.6444+0.107276*t-1.3839e-3*t**2)         &
           +s**(3./2.)*(0.17709-4.0772e-3*t+5.3539e-5*t**2)

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                      TEMP_PROFILE
!   computes diffusive mixing of temperature profile
!   this is the implicit crank-nicholson version -- sept 5 1996
!   implimented for speed and to minimize possible CFL problems
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine temp_profile( iwater, qbot, qw, t, sw, lnet, Qe, &
                               Qh, de, depth, salty)
      implicit none
      include 'Tanganyika.inc'
      real qbot,qw,sw,lnet,Qe,Qh,area_2,t1,cnextra,swtop,    &
           dist12,area_1,top,bot
      integer az,k,iwater,ktop

      parameter (az = 1000 ) ! length of arrays first dimensioned here
      real t(depth,2), de(depth), dnsty(depth)
      real cpz(az), z(az), zhalf(az) ! spec heat, dz, and dz.5
      real a(az), b(az), c(az), d(az)  ! arrays for tridia matrix
      real told(az), tnew(az)
      real salty(max_dep,n_trace)

! z is an array of delz with depth
! zhalf is array of delz between i and i+1
! so, zhalf(3) = dist from mid of 3 to mid of 4
! zhalf = .5*(dz(i)+dz(i+1))
! here all 1 except for zhalf(1) = .8

      dist12=surf*0.5+dz*0.5 ! dist for grad t calc between
      ktop = max_dep-depth  ! index value of lake surface

      do k= 1,depth ! init density, spec_ht, zhalf through col
        t(k,2)=t(k,1)
        told(k)=t(k,1)
        call density(t(k,1),salty(k,n_trace),dnsty(k))
        call specheat(t(k,1),salty(k,n_trace),cpz(k))
        z(k)=dz
        zhalf(k)=dz
      enddo

        z(1)=surf                ! make top slice z = surf
        zhalf(1)=0.5*(z(1)+z(2)) ! calc zhalf between 1 and 2

!-----------------------------------------------------------------------
!   2. calculate d array, which includes heating terms
!      so seperate for ice and water fraction of cell
!     cnextra is extra term in d array for crank-nicholson
!------------------------------------------------------------------------

!   2.1   First do top slice of column
           k = 1
           area_1 =(area(k+ktop) + area(k+ktop)) /2.
! 0.6m above base of dz is almost = to where basin area is assigned
           area_2 =(area(k+ktop) + area(k+ktop+1)) / 2.
          if (iwater.eq.1) then  ! open water calculation
           t1 = sw*beta +   &
                (1.-beta)*sw*(1.-exp(-eta*surf))*area_2/area(k+ktop) +  &
                (Lnet+Qe+Qh) * area_1/area(k+ktop)
          else                   ! calculation beneath ice
           t1 = qbot*beta +       &
                (1.-beta)*qbot*(1.-exp(-eta*surf))*area_2/area(k+ktop) - &
                qw * area_1/area(k+ktop)
          endif

        cnextra = 0.5 * area_2/area(k+ktop) * & ! for k = 1, top slice in column
       ((de(1)/zhalf(1)) * (dt/z(1)) * (t(2,1) - t(1,1)))

        d(1) = t(1,1)+t1*dt/((1.e3+dnsty(1))*cpz(1)*z(1))+cnextra

!   2.2   Remainder of water column, excluding top and bottom slice
       swtop = qbot   ! use ice value
       if (iwater.eq.1) swtop = sw ! set to water value if wat frac.
       do  k=2,depth-1           ! loop through water column
         top = (surf+(k-2)*dz)
         bot = (surf+(k-1)*dz)
         area_1 =(area(k+ktop-1) + area(k+ktop)) / 2.
         area_2 =(area(k+ktop) + area(k+ktop+1)) / 2.
         t1 = (1.-beta)*swtop*  &
              ((area_1*exp(-eta*top)-area_2*exp(-eta*bot))/area(k+ktop))

!      terms below reversed for C-N

         cnextra = 0.5 *1./area(k+ktop)*  &
             (((de(k)/zhalf(k))*(dt/z(k))*(t(k+1,1)-t(k,1))) *area_2 -  &
             ((de(k-1)/zhalf(k-1))*(dt/z(k))*(t(k,1)-t(k-1,1)))*area_1)
         d(k) = t(k,1)+t1*dt/((1.e3+dnsty(k))*cpz(k)*z(k))+cnextra

       enddo                   ! do next slice in column

!    2.3  Bottom slice of the water column
!         No mud heating
      k=depth
      area_1 =(area(k+ktop-1) + area(k+ktop)) / 2.
      swtop = qbot   ! use ice value
      if (iwater.eq.1) swtop = sw ! set to water value if wat frac.
      top = (surf+(k-2)*dz)
      t1 = (1.-beta)*swtop*(exp(-eta*top))*area_1/area(k+ktop)
      cnextra = 0.5 * area_1/area(k+ktop) * &  ! for k = depth
       ((de(k-1)/zhalf(k-1)) * (dt/z(k)) * (t(k,1) - t(k-1,1)))
      d(k) = t(k,1)+t1*dt/((1.e3+dnsty(k))*cpz(k)*z(k))+cnextra

!-----------------------------------------------------------------------
!   3. calculate arrays for tridiagonal matrix, top and bot slice seperate
!------------------------------------------------------------------------

      k = 1  ! do top slice
      area_2 =(area(k+ktop) + area(k+ktop+1)) / 2.
      b(1) = -0.5 * (de(1)/zhalf(1)) *   &
              (dt/z(1)) * area_2/area(k+ktop)
      a(1) = 1. - b(1)  ! no factor of 0.5 here, already done

      do  k = 2,depth-1 ! do all but top and bottom slices
        area_1 =(area(k+ktop-1) + area(k+ktop)) / 2.
        area_2 =(area(k+ktop) + area(k+ktop+1)) / 2.
        b(k) = -0.5 * (de(k)/zhalf(k)) *  &
              (dt/z(k))*area_2/area(k+ktop)
        c(k) = -0.5 * (de(k-1)/zhalf(k-1)) *  &
              (dt/z(k))*area_1/area(k+ktop)
        a(k) = 1. - b(k) - c(k) ! no 0.5 here, already done
      enddo
      k=depth  ! do bottom slice
      area_1 =(area(k+ktop-1) + area(k+ktop)) / 2.
      c(depth) = -0.5 * (de(depth)/zhalf(depth) ) *  &
              (dt/z(depth)) * area_1/area(k+ktop)
      a(depth) = 1. - c(depth)

!-------------------------------------------------------------
!    4. call matrix solution, and reset temp and density arrays
!-------------------------------------------------------------

       call tridiag_solve(1,1,depth,c,a,b,d,tnew)
       do k = 1, depth  ! change temps and densities
          t(k,2) = tnew(k)
          t(k,1) = t(k,2)
          call density (t(k,2),salty(k,n_trace),dnsty(k))
       enddo

       return
       end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                   TENDENCY
!     computes/interpolates the tendency between input reads
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !   call tendency (j,nsteps,ta_in,ta_i,qa_in,qa_i,ua_in,ua_i,   &
      !                  sw_in,sw_i,rlwd_in,rlwd_i,ps_in,ps_i,        &
      !                  rh_in,rh_i,prec_in,prec_i,deutprec_in, &
      !                  deutprec_i,o18prec_in, o18prec_i, &
      !                  runin_in,runin_i,deutrun_in, deutrun_i,&
      !                  o18run_in,o18run_i)
      !Ashling: added isotopes to tendency (comment above shows call)

      subroutine tendency(li,idtend,ta_in,ta_i,qa_in,qa_i,  &
                          ua_in,ua_i,sw_in,sw_i,rlwd_in,rlwd_i,  &
                          ps_in,ps_i,rh_in,rh_i,prec_in,prec_i,  &
                          deutprec_in, deutprec_i, o18prec_in,&
                          o18prec_i,runin_in,runin_i, deutrun_in,&
                          deutrun_i,o18run_in,o18run_i)
      implicit none
      include 'Tanganyika.inc'
      real ta_in,ta_i,qa_in,qa_i,ua_in,ua_i,sw_in,sw_i,rh_in,rh_i, &
           rlwd_in,rlwd_i,ps_in,ps_i,prec_in,prec_i,runin_in,runin_i, &
           deutrun_in, deutprec_in, o18prec_in, o18run_in, &
           deutrun_i, deutprec_i, o18prec_i, o18run_i !Ashling
      dimension ta_in(2),qa_in(2),ua_in(2),sw_in(2),rh_in(2),  &
                rlwd_in(2),ps_in(2),prec_in(2),runin_in(2), &
                deutrun_in(2), deutprec_in(2), o18prec_in(2), &
                o18run_in(2) !Ashling
      integer li,idtend

      ta_i = ta_in(1)+(li-1)*(ta_in(2)-ta_in(1))/real(idtend)
      qa_i = qa_in(1)+(li-1)*(qa_in(2)-qa_in(1))/real(idtend)
      ua_i = ua_in(1)+(li-1)*(ua_in(2)-ua_in(1))/real(idtend)
      sw_i = sw_in(1)+(li-1)*(sw_in(2)-sw_in(1))/real(idtend)
      rlwd_i = rlwd_in(1)+(li-1)*(rlwd_in(2)-rlwd_in(1))/real(idtend)
      ps_i = ps_in(1)+(li-1)*(ps_in(2)-ps_in(1))/real(idtend)
      rh_i = rh_in(1)+(li-1)*(rh_in(2)-rh_in(1))/real(idtend)

      !Begin Ashling
      if(deutflag) then
      deutprec_i=deutprec_in(1)+(li-1)* &
                 (deutprec_in(2)-deutprec_in(1))/real(idtend)
      deutrun_i=deutrun_in(1)+(li-1)* &
                 (deutrun_in(2)-deutrun_in(1))/real(idtend)
      end if
      if(o18flag) then
      o18prec_i=o18prec_in(1)+(li-1)* &
                 (o18prec_in(2)-o18prec_in(1))/real(idtend)
      o18run_i=o18run_in(1)+(li-1)* &
                 (o18run_in(2)-o18run_in(1))/real(idtend)
      end if
      !End Ashling

      prec_i = prec_in(1)/real(idtend)  ! accumulations, not rates
      runin_i = runin_in(1)/real(idtend)

      if (li.eq.idtend) then
         ta_in(1)=ta_in(2)
         qa_in(1)=qa_in(2)
         ua_in(1)=ua_in(2)
         sw_in(1)=sw_in(2)
         rlwd_in(1)=rlwd_in(2)
         ps_in(1)=ps_in(2)
         rh_in(1)=rh_in(2)
         prec_in(1)=prec_in(2)
         runin_in(1)=runin_in(2)
         !begin Ashling
         deutprec_in(1)=deutprec_in(2)
         deutrun_in(1)=deutrun_in(2)
         o18prec_in(1)=o18prec_in(2)
         o18run_in(1)=o18run_in(2)
         !end Ashling
      endif

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                   TRACER_AVG
!     average water and ice for fractional cover
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine tracer_avg (depth, fracprv)

      implicit none
      include 'Tanganyika.inc'
      real fracprv
      integer j,i_tr

      do j=1,depth
         do i_tr = 1,num_tra
            trace (j,i_tr) = (1.-fracprv)*trace(j,i_tr) + &
                                 fracprv*trace_i(j,i_tr)
         enddo
      enddo

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                TRACER_MIXER
!  mixes lake layers due to density instabilities
!  temperature as well as tracers(isotopes, salinity) get mixed
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine tracer_mixer (t, dnsty, depth, salty, mixdep, iwater)

      implicit none
      include 'Tanganyika.inc'
      integer k,mixprev,i_tr,iwater,m,mixdep,kk,k2
      real avet,avev,avev_tr,cp,vol,vol_tr,tav,densnew,rho_max
      real t(depth,2), dnsty(depth), tr_work(max_dep,n_trace)
      real tr_av (n_trace), ave_tr (n_trace), salty(max_dep,n_trace)

      do k=1,depth
         call density(t(k,1),salty(k,n_trace),dnsty(k))
      enddo

      mixprev = 1 ! top depth of local instability, may not be 1

      do  k= 1,depth
        do  i_tr = 1, num_tra
          if (iwater.eq.1) tr_work(k,i_tr) = trace (k,i_tr)  ! open water calc
          if (iwater.eq.0) tr_work(k,i_tr) = trace_i (k,i_tr)  ! ice calc
        enddo
      enddo

!-----------------------------------------------------------------------
!   1. check for density instability at each slice in water column
!------------------------------------------------------------------------
 9     continue            ! if a new instability created by mixer
       do  k=1,depth-1
        avet=0.0
        avev=0.0

        do i_tr = 1, num_tra
          ave_tr (i_tr) = 0.
        enddo
        avev_tr = 0.    ! out of tracer loop, same for all tracers

      if (dnsty(k).gt.dnsty(k+1)) then

!  if density in layer above greater than density in layer below, keep
!  track of the depth to which mixing should occur (mixdep)

        if (mixprev.eq.1.0 .and.(k+1).gt.mixdep) &  ! reassigns mix
            mixdep=k+1 ! if top mix and dep > wat dep for ice

!-----------------------------------------------------------------------
!   2.  sum and average temp,tracer,volume,density from mixprev to k+1
!------------------------------------------------------------------------

        do m = mixprev, k+1 ! mix from top of instab to bot of inst
           call specheat (t(m,2),salty(m,n_trace),cp)
           if (m.eq.1) then
              vol = surf*(1.e3+dnsty(m))*cp*area(m+max_dep-depth) ! joules / deg
              vol_tr = surf * area(m+max_dep-depth)
!  calculate the total heat in each layer of unstable water (vol) and the
!  total volume of each layer (vol_tr)
           else
              vol = dz*(1.e3+dnsty(m))*cp*area(m+max_dep-depth)   ! joules / deg
              vol_tr = dz * area(m+max_dep-depth)
!  same as above, but using dz instead of surf to reflect the fact that we
!  are in a layer below the surface
           endif
           avet=avet+t(m,2)*vol
           avev=avev+vol

!  sum up the total amount of heat contained within the slab of unstable water
!  by adding all the slices together, then sum up the total volume of the
!  slab and the total amount of any tracers (see below)

           do   i_tr = 1, num_tra
             ave_tr(i_tr)=ave_tr(i_tr)+tr_work(m,i_tr)*vol_tr
           enddo
           avev_tr = avev_tr + vol_tr   ! sum of volume

        enddo  ! end summing

!  calculate the new average temperature and tracer concentrations for the
!  mixed layer

        Tav=avet/avev
        do i_tr = 1,num_tra
          tr_av (i_tr) = ave_tr(i_tr) / avev_tr
        enddo

!  next, given the new average temperature and salinity, calculate the new
!  density profile for the water
        if (s_flag) then    ! need to use mixed value, which is the last tracer
          call density (Tav,tr_av(num_tra),densnew)
        else    ! salinity not changing, just use surface value
          call density (Tav,salty(1,n_trace),densnew)
        endif

!   2.1 check to make this has not generated new instabilities
!       above mixprev

        rho_max = -50.0   !

!	write(*,*)'mixprev is',mixprev

        do kk=1,mixprev-1   ! find maximum density above mixed section
         if (dnsty(kk).gt.rho_max ) rho_max=dnsty(kk)
        enddo

!-----------------------------------------------------------------------
!   3.  adjust temp,tracers and density in the mixed part of column
!------------------------------------------------------------------------

        do k2 = mixprev, k+1
           t(k2,2)=Tav
           do i_tr = 1,num_tra
             tr_work(k2,i_tr) = tr_av (i_tr)
           enddo
           dnsty(k2)=densnew
        enddo

!   3.1 compare rho_max to adjusted new densities
        if (rho_max.gt.(densnew)) then   ! need to go back to top
         mixprev = 1                   ! reset to top of column
         go to 9                       ! start checking column again
        endif                          ! if new instability created

       else         ! if density at K not > k + 1
        mixprev=k+1
       endif        ! if density at k > k+1

      enddo  !   end looping through water column searching for inst.

      do k = 1, depth
             t(k,1)=t(k,2)

        do i_tr = 1,num_tra
          if (iwater.eq.1) trace(k,i_tr)=tr_work(k,i_tr) ! open water calc
          if (iwater.eq.0) trace_i(k,i_tr)=tr_work(k,i_tr) ! ice calc
      enddo

!  once mixing is completely done, calculate the new density profile for the
!  mixed water

        call density (t(k,1),salty(k,n_trace),dnsty(k))
      enddo

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               TRACER_PROFILE
!     computes diffusive mixing of passive tracers
!     this is a crank-nicholson version -- sept 5 1996
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine tracer_profile (de, depth, iwater)

      implicit none
      include 'Tanganyika.inc'
      real dist12,area_1,area_2,cnextra
      integer az,i_tr,k,i_water,iwater,ktop
      parameter (az = 1000) ! length of arrays first dimensioned here

      real de(max_dep)
      real z(az), zhalf(az) ! dz, and dz.5
      real a(az), b(az), c(az), d(az)  ! arrays for tridia matrix
      real tnew(az)
      real tr_work (max_dep)

! z is an array of delz with depth
! zhalf is array of delz between i and i+1
! so, zhalf(3) = dist from mid of 3 to mid of 4
! zhalf = .5*(dz(i)+dz(i+1))
! here all 1 except for zhalf(1) = .8

!-----------------------------------------------------------------------
!   0.  Loop through the different tracers
!------------------------------------------------------------------------

      do 1 i_tr = 1, num_tra

!-----------------------------------------------------------------------
!   1.  Place correct tracer into work arrray
!------------------------------------------------------------------------

      dist12=surf*0.5+dz*0.5 ! dist for grad t calc between
      ktop = max_dep-depth

      do 10 k= 1,depth ! init density, spec_ht, zhalf through col
        z(k)=dz
        zhalf(k)=dz
        if (iwater.eq.1) tr_work(k) = trace (k,i_tr)  ! open water calc
        if (iwater.eq.0) tr_work(k) = trace_i (k,i_tr)  ! ice calc
 10   continue

      z(1)=surf                ! make top slice z = surf
      zhalf(1)=0.5*(z(1)+z(2)) ! calc zhalf between 1 and 2

!-----------------------------------------------------------------------
!   2. calculate d array, which includes heating terms
!      so seperate for ice and water fraction of cell
!     cnextra is extra term in d array for crank-nicholson
!------------------------------------------------------------------------

      k = 1   ! top slice
      area_1 =(area(k+ktop) + area(k+ktop)) / 2.
      area_2 =(area(k+ktop) + area(k+ktop+1)) / 2.

      cnextra = 0.5 * area_2/area(k+ktop) * & ! for k = 1, top slice in column
       ((de(1)/zhalf(1)) * (dt/z(1))*(tr_work(2) - tr_work(1)))

      d(1) = tr_work(1)+cnextra

      do 11 k=2,depth-1           ! loop through water column

        area_1 =(area(k+ktop-1) + area(k+ktop)) / 2.
        area_2 =(area(k+ktop)  + area(k+ktop+1)) / 2.

        cnextra = 0.5 *1./area(k+ktop)*                     &
                 (((de(k)/zhalf(k))*(dt/z(k))*         &
                 (tr_work(k+1)-tr_work(k)))*area_2 -   &
                 ((de(k-1)/zhalf(k-1))*(dt/z(k))*      &
                 (tr_work(k)-tr_work(k-1)))*area_1)

        d(k) = tr_work(k)+cnextra

 11   continue

       k = depth  !  bottom slice
       area_1 =(area(k+ktop-1) + area(k+ktop)) / 2.

       cnextra = 0.5 * area_1/area(k+ktop) *   & ! for k = depth
       ((de(k-1)/zhalf(k-1))*(dt/z(k))*(tr_work(k)-tr_work(k-1)))

       d(k) = tr_work(k)+cnextra

!-----------------------------------------------------------------------
!   3. calculate arrays for tridiagonal matrix, top and bot slice seperate
!------------------------------------------------------------------------

      k = 1  ! top slice
      area_2 =(area(k+ktop) + area(k+ktop+1)) / 2.

      b(1) = -0.5 * (de(1)/zhalf(1)) *   &
              (dt/z(1)) * area_2/area(k+ktop)
      a(1) = 1. - b(1)  ! no factor of 0.5 here, already done

      do 13 k = 2,depth-1  ! all but top and bottom slices

         area_1 =(area(k+ktop-1) + area(k+ktop)) / 2.
         area_2 =(area(k+ktop) + area(k+ktop+1)) / 2.

         b(k) = -0.5 * (de(k)/zhalf(k)) *   &
              (dt/z(k))*area_2/area(k+ktop)
         c(k) = -0.5 * (de(k-1)/zhalf(k-1)) *  &
              (dt/z(k))*area_1/area(k+ktop)
         a(k) = 1. - b(k) - c(k) ! no 0.5 here, already done

 13   continue

      k=depth  ! bottom slice
      area_1 =(area(k+ktop-1) + area(k+ktop)) / 2.

      c(depth) = -0.5 * (de(depth)/zhalf(depth)) *   &
              (dt/z(depth)) * area_1/area(k+ktop)
      a(depth) = 1. - c(depth)

!-------------------------------------------------------------
!    4. call matrix solution, and reset temp and density arrays
!-------------------------------------------------------------

      call tridiag_solve(1,1,depth,c,a,b,d,tnew)

      do 40 k = 1, depth  ! change temps and densities
          tr_work(k) = tnew(k)
        if (iwater.eq.1) trace(k,i_tr) = tr_work(k) ! open water calc
        if (iwater.eq.0) trace_i(k,i_tr) = tr_work(k) ! ice calc
 40   continue


!-------------------------------------------------------------
 1    continue    ! go back and do the next tracer

      return
     end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                    TRIDIAG_SOLVE
!     to compute the solution of many tridiagonal linear systems.
!     used to solve system of equations for eddy diffusion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine tridiag_solve (ns, nd, ne, a, b, c, y, x)

!      ns ..... the number of systems to be solved.
!      nd ..... first dimension of arrays (ge ns).
!      ne ..... the number of unknowns in each system.
!               this must be > 2. second dimension of arrays.
!      a ...... the subdiagonals of the matrices are stored
!               in locations a(j,2) through a(j,ne).
!      b ...... the main diagonals of the matrices are stored
!               in locations b(j,1) through b(j,ne).
!      c ...... the super-diagonals of the matrices are stored in
!               locations c(j,1) through c(j,ne-1).
!      y ...... the right hand side of the equations is stored in
!               y(j,1) through y(j,ne).
!      x ...... the solutions of the systems are returned in
!               locations x(j,1) through x(j,ne).
!      alpha .. work array dimensioned alpha(nd,ne)
!      gamma .. work array dimensioned gamma(nd,ne)
!      history: based on a streamlined version of the old ncar
!               ulib subr trdi used in the phoenix climate
!               model of schneider and thompson (j.g.r., 1981).
!               revised by starley thompson to solve multiple
!               systems and vectorize well on the cray-1.
!               later revised to include a parameter statement
!               to define loop limits and thus enable cray short
!               vector loops.
!      algorithm: lu decomposition followed by solution.
!                 note: this subr executes satisfactorily
!                 if the input matrix is diagonally dominant
!                 and non-singular.  the diagonal elements are
!                 used to pivot, and no tests are made to determine
!                 singularity. if a singular or numerically singular
!                 matrix is used as input a divide by zero or
!                 floating point overflow will result.

!     last revision date:      4 february 1988

      implicit none
      integer ns,ne,nd,nm1,j,i,ib
      real a(nd,ne), b(nd,ne), c(nd,ne), y(nd,ne),  &
           x(nd,ne), alpha(nd,ne), gamma(nd,ne)

      nm1 = ne-1

!        obtain the lu decompositions.

      do j=1,ns
         alpha(j,1) = 1./b(j,1)
         gamma(j,1) = c(j,1)*alpha(j,1)
      enddo
      do i=2,nm1
         do j=1,ns
            alpha(j,i) = 1./(b(j,i)-a(j,i)*gamma(j,i-1))
            gamma(j,i) = c(j,i)*alpha(j,i)
         enddo
      enddo

      do j=1,ns
         x(j,1) = y(j,1)*alpha(j,1)
      enddo
      do i=2,nm1
         do j=1,ns
            x(j,i) = (y(j,i)-a(j,i)*x(j,i-1))*alpha(j,i)
         enddo
      enddo
      do j=1,ns
         x(j,ne) = (y(j,ne)-a(j,ne)*x(j,nm1))/   &
                   (b(j,ne)-a(j,ne)*gamma(j,nm1))
      enddo
      do i=1,nm1
         ib = ne-i
         do j=1,ns
            x(j,ib) = x(j,ib)-gamma(j,ib)*x(j,ib+1)
         enddo
      enddo

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  WATER_BALANCE
!    calculates the water balance of the lake
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine water_balance (depth,d_frac,prec,evap1,evap2,snowmelt,&
                                runin, runout, surf_a, n_slice, isave_d)
      implicit none
      include 'Tanganyika.inc'
      real prec, runin, run_len, surf_a, elevsurf, run_s,   &
           d_level, remain, d_area, s_surf, r_salt, runout, &
           evap1, evap2, snowmelt
      integer isave_d, i_surf, n_slice, k_tr, i_tr

      run_s = 0.0       ! runoff salinity

!==============================================================
!      1. convert runoff input/output volume to rate
!==============================================================

      i_surf = max_dep - depth + 1 ! lk_area index of water surface
      if (i_surf.eq.1) then
         d_area = 0.0              ! lake at sill
      else
         d_area = area(i_surf-1) - area(i_surf)
      end if
      surf_a = d_frac * (d_area/dz) + area(i_surf) ! actual surf area
      run_len = (runin*(b_area-surf_a)/surf_a)

!==============================================================
!      2. calculate change in lake level
!==============================================================

       if (wb_flag) then
          d_level = (run_len+prec-evap1/1.e3*dt-evap2/1.e3*dt) -  &
                    snowmelt*rhosnow/rhowat
          d_frac  = d_frac + d_level
          remain = mod(d_frac,dz)          !remainder over whole dzs
          n_slice = int(d_frac - remain )  !number of new slices
          d_frac = remain
          depth = depth + n_slice           !new depth in integer slices

          if (depth .ge. max_dep.and.d_frac.gt.0.) then
            runout = d_frac
            d_frac = 0.0
            depth=max_dep
          else
            runout = 0.0
          endif
       endif

!==============================================================
!      2.5 calculate surface salt balance
!==============================================================

       if (s_flag) then
          s_surf = trace(1,n_trace) * surf  ! volume of salt (ppt*m)
          r_salt = run_len*run_s*dt/1.e3    ! volume of runoff salt (ppt*m)
          s_surf = (s_surf + r_salt)*( 1./(surf + d_level))
          trace(1,n_trace) = s_surf
       endif

!==============================================================
!      3. adjust change in lake level, add new layers
!==============================================================

       if (wb_flag) then

       elevsurf=(basedep+depth+d_frac)
       if (depth.ne.isave_d) then  ! adjust area, depth has changed
!        print*,'adjusting depth new elevation ',elevsurf,' m'
         if (depth.lt.0.) stop ' Lake is dry '
       endif

       if (depth .ne. isave_d) then
         do  k_tr = 1, n_slice
           do  i_tr = 1, num_tra
            trace(isave_d+k_tr,i_tr) = trace (isave_d,i_tr)
           enddo
         enddo
       endif

       endif

      return
      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 ZERO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine zero (t2w,t2i,q2w,q2i,u2w,u2i,evapw,evapi,qhw,qhi, &
                        qew,qei,lnetw,lneti,luw,lui,sww,swi)
      implicit none
      real t2w,t2i,q2w,q2i,u2w,u2i,evapw,evapi,qhw,qhi,             &
           qew,qei,lnetw,lneti,luw,lui,sww,swi

      t2w=0.0
      t2i=0.0
      q2w=0.0
      q2i=0.0
      u2w=0.0
      u2i=0.0
      evapw=0.0
      evapi=0.0
      qhw=0.0
      qhi=0.0
      qew=0.0
      qei=0.0
      lnetw=0.0
      lneti=0.0
      luw=0.0
      lui=0.0
      sww=0.0
      swi=0.0

      return
      end
