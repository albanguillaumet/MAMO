# This file contains some data associated with the program 

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# I - DAILY RATES OF MALARIA TRANSMISSION RISK

# Our data consist of the daily rate of malaria transmission risk at 11 stations sampled at 100 meter intervals between 2000 and 1000 m ( i.e, grad = c(2000, 1000) ). 

# 1 - Year divided into three seasons (Apapane and Iiwi)

# 1.1 - Variables corresponding to grad = c(2000, 1000)

# The current daily rate for the IIWI and APAP breeding season (September-April)
alpha.b_nc.11 = c(0.00000e+00, 0.00000e+00, 0.00000e+00, 3.21323e-06, 6.42646e-06, 1.40851e-03, 3.37271e-03, 5.89904e-03, 8.98748e-03, 1.29177e-02, 1.34780e-02)
# The current daily rate for the IIWI and APAP first period of the non-breeding season (May-June)
alpha.nb.1_nc.11 = c(0.00000e+00, 0.00000e+00, 0.00000e+00, 1.61291e-08, 3.22581e-08, 7.02270e-04, 1.68541e-03, 2.94945e-03, 4.49439e-03, 6.46068e-03, 6.74157e-03)
# The current daily rate for the IIWI and APAP second period of the non-breeding season (July-August)
alpha.nb.2_nc.11 = c(0.00000e+00, 0.00000e+00, 0.00000e+00, 6.92295e-05, 1.38459e-04, 2.92186e-03, 6.84632e-03, 1.19118e-02, 1.81184e-02, 2.60054e-02, 2.71084e-02 )
# The daily rate in 2100 for the IIWI and APAP breeding season (September-April)
alpha.b.2100_nc.11 = c(6.426460e-06, 3.467123e-04, 6.869982e-04, 1.027284e-03, 1.367570e-03, 4.906265e-03, 1.013395e-02, 1.705063e-02, 2.565630e-02, 3.648775e-02, 3.780066e-02)
# The daily rate in 2100 for the IIWI and APAP first period of the non-breeding season (May-June)
alpha.nb.1.2100_nc.11 = c(3.225810e-08, 6.282144e-05, 1.256106e-04, 1.883998e-04, 2.511890e-04, 4.162655e-03, 9.688946e-03, 1.683006e-02, 2.558600e-02, 3.670766e-02, 3.825341e-02)
# The daily rate in 2100 for the IIWI and APAP second period of the non-breeding season (July-August)
alpha.nb.2.2100_nc.11 = c(0.0001384590, 0.0005414685, 0.0009444780, 0.0013474875, 0.0017504970, 0.0109072012, 0.0240766864, 0.0412589527, 0.0624540000, 0.0892743571, 0.0928057515)

# 1.2 - The same variables corresponding to grad = c(1900, 1000), i.e. 10 stations sampled at 100 meter intervals between 1900 and 1000 m
alpha.b = alpha.b_nc.11[-1] 
alpha.nb.1 = alpha.nb.1_nc.11[-1] 
alpha.nb.2 = alpha.nb.2_nc.11[-1] 
alpha.b.2100 = alpha.b.2100_nc.11[-1] 
alpha.nb.1.2100 =  alpha.nb.1.2100_nc.11[-1]
alpha.nb.2.2100 = alpha.nb.2.2100_nc.11[-1] 

# 2 - A single season (other species)

# 2.1 - Variables corresponding to grad = c(2000, 1000)

# The current daily rate
alpha.1_nc.11 = c(0.000000e+00, 0.000000e+00, 0.000000e+00, 1.351005e-05, 2.702011e-05, 1.538971e-03, 3.661105e-03,  6.393424e-03, 9.735927e-03, 1.398763e-02, 1.459038e-02)
# The daily rate in 2100
alpha.1.2100_nc.11 = c(2.702011e-05, 3.298001e-04, 6.325801e-04, 9.353600e-04, 1.238140e-03, 5.823982e-03, 1.249179e-02, 2.124156e-02, 3.207330e-02, 4.574940e-02, 4.749088e-02)

# 2.2 - The same variables corresponding to grad = c(1900, 1000)

alpha.1 = alpha.1_nc.11[-1] 
alpha.1.2100 = alpha.1.2100_nc.11[-1]

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# II - DISTRIBUTION OF RESOURCES (OHIA BLOOM) DURING THE NON-BREEDING SEASON (FOR IIWI AND APAPANE)

# Our data consist of a blooming index for ohia (Metrosideros polymorpha) trees at 10 stations sampled at 100 meter intervals between 1900 and 1000 m ( i.e, grad = c(1900, 1000) ). 
# see Kuntz (2008) for original data source, and Guillaumet et al. (submitted) for additional inferences.

# 1 - Data corresponding to grad = c(1900, 1000)

# The bloom index during the first period of the non-breeding season (May-June) in 2003
K.nb.1.2003 = c(2.45556,  3.23333,  4.01111,  4.78889,  5.56667,  6.05000,  9.81667, 17.52500, 19.12000, 17.33333)
# The bloom index during the first period of the non-breeding season (May-June) in 2004
K.nb.1.2004 = c(5.51333, 9.58000, 13.64667, 17.71333, 21.78000, 10.00000,  7.86667,  4.60000, 10.14000, 11.78571)
# The average bloom index during the first period of the non-breeding season (May-June) in 2003 and 2004
K.nb.1.avg = 0.5 * (K.nb.1.2003 + K.nb.1.2004)
# The bloom index during the second period of the non-breeding season (July-August; average of 2003 and 2004) 
K.nb.2 = c(4.06815, 7.07556, 10.08296, 13.09037, 16.09778, 19.16667,  7.47778, 4.55000, 4.16667,  2.36667)

# 2 - We also create the following variables corresponding to grad = c(2000, 1000)
# We assume that the value at 2000 m equals the value at 1900 m following reforestation. 

K.nb.1.2003_nc.11 = c(K.nb.1.2003[1], K.nb.1.2003)
K.nb.1.2004_nc.11 = c(K.nb.1.2004[1], K.nb.1.2004)
K.nb.1.avg_nc.11 = c(K.nb.1.avg[1], K.nb.1.avg)
K.nb.2_nc.11 = c(K.nb.2[1], K.nb.2)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# III - SPECIES-SPECIFIC DENSITY ( # pairs / ha)  ALONG THE ELEVATIONAL GRADIENT grad = c(1900, 1000)

# 1 - IIWI in Hamakua

# estimates obtained based on DISTANCE-corrected estimates; see Guillaumet et al. (submitted) for details

y.obs.IIWI.1 = c(5.39163173, 4.84949943, 5.11686774, 4.70852507, 3.66411526, 1.75865820, 0.24413093, 0.00000000, 0.00887159, 0.00000000)

# 2 - AT HAKALAU 

# The data was obtained by combining Hamakua data at 1900-1500 m and Kau data at 1400-1000 m, while ensuring that the maximum density was the one observed at Hakalau NWR (see Guillaumet et al. in prep. for details)
y.obs.AKIP = c(0.057857604, 0.080000000, 0.074552759, 0.044400641, 0.043278985, 0.035770454, 0.000000000, 0.009895678, 0.000000000, 0.000000000)
y.obs.AKEP = c(0.97000000, 0.65407215, 0.71964077, 0.62545589, 0.32342036, 0.04394132, 0.17235153, 0.21069672, 0.30465185, 0.00000000)
y.obs.HCRE = c(0.80292661, 0.68062835, 0.92477397, 1.17000000, 0.95131759, 0.26290216, 0.32538981, 0.39778332, 0.07189565, 0.00000000)
y.obs.IIWI = c(8.250000, 7.463660, 7.540451, 7.428597, 6.117269, 3.039700, 0.000000, 0.000000, 0.000000, 0.000000)
y.obs.ELEP = c(1.0168305, 1.4506461, 1.7200000, 1.4398433, 1.6909516, 1.3888620, 1.2987445, 1.7200000, 1.2195829, 0.8678246)
y.obs.OMAO = c(0.8815179, 0.9194332, 0.8806342, 0.9078896, 0.9200000, 0.8914893, 0.7991310, 0.7399213, 0.6638720, 0.6659613)
y.obs.APAP = c(4.8629363, 4.9600000, 4.0586814, 4.2813933, 4.4315999, 4.9357569, 2.1841874, 1.0456216, 0.2788324, 0.0000000)
y.obs.HAAM = c(5.207966, 4.406453, 5.240000, 4.302408, 4.066120, 3.526126, 4.359680, 2.043600, 3.269760, 3.633067)
