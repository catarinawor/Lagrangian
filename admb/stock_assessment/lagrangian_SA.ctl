## ———————————————————————————————————————————————————————————————————————————————————— ##
## CONTROL FILE - stolen from iSCAM                                                     ##
## ———————————————————————————————————————————————————————————————————————————————————— ##
## CONTROLS FOR LEADING PARAMETERS                                                      ##
##  Prior descriptions:                                                                 ##
##                      -0 uniform      (0,0)                                           ##
##                      -1 normal       (p1=mu,p2=sig)                                  ##
##                      -2 lognormal    (p1=log(mu),p2=sig)                             ##
##                      -3 beta         (p1=alpha,p2=beta)                              ##
##                      -4 gamma        (p1=alpha,p2=beta)                              ##
## ———————————————————————————————————————————————————————————————————————————————————— ##
# npar 
8
##
## ival         lb      ub      phz     prior   p1      p2        #parameter            ##
## ———————————————————————————————————————————————————————————————————————————————————— ##
1.09861	0	1.8	3	0	0	1.8	#log_mo
-2.3	-3	-0.01	3	0	-3	-0.01	#log_cvPos
1.6	1	2.1	3	2	0.47	0.05	#log_maxPos50
0.5	-0.4	1.4	3	0	-0.4	1.4	#log_maxPossd
0.980579	-0.7	3	1	1	1	0.2	#log_Ro 
0.862	0.2	1	1	3	9.909	2.959	#h
2	-0.7	6	2	0	-0.7	6	#log_avgrec
10.491	0.01	10	1	2	1.6	0.5	#sigma_r 
## ———————————————————————————————————————————————————————————————————————————————————— ##
##
## wt_ival
 -1.6815 -1.6722 1.6044 0.5671 -1.7991 1.4167 -0.6429 -1.7381 1.1481 0.968 0.2378 0.482 0.0969 0.6133 2.472 -1.0325 0.0921 -2.7704 0.2381 -2.4648 0.7649 0.5756 -2.8746 1.7223 0.1782 2.7192 -0.6564 0.3729 -0.1946 2.4267	 -1.6815 -1.6722 1.6044 0.5671 -1.7991 1.4167 -0.6429 -1.7381 1.1481 0.968 0.2378 0.482 0.0969 0.6133 2.472 -1.0325 0.0921 -2.7704 0.2381 -2.4648
##
