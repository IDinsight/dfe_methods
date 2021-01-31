
* Specify output folder
global WD "C:/Users/tofis/Dropbox (IDinsight)/IE Design for DFEs/Code/03_Output"


* Standard Approach: 0.25 sd effect size, mean 0.2, two-sided, alpha .05 
power twomeans 0.2, nratio(1) /// 
	power(0.8) alpha(0.05) sd(1) diff(0.05)

* Test 1: 0.25 sd effect size, mean 0, onesided, alpha .05 
power twomeans 0, onesided nratio(1) /// 
	power(0.8) alpha(0.05) sd(1) diff(0.25)
	

* Test 2: 0.05 sd effect size, mean of .2, onesided, alpha 0.2 

power twomeans 0.2, onesided nratio(1) /// 
	power(0.8) alpha(0.2) sd(1) diff(0.05)	
	
	
	
* Conservative Test	

power twomeans 0.2, onesided nratio(1) /// 
	power(0.8) alpha(0.05) sd(1) diff(0.05)	
	
power twomeans 0.2, onesided nratio(1) /// 
	power(0.8) alpha(0.2) sd(1) diff(0.05)	

	
* Figure 2 - Trade-offs b/w sample size and p-values

power twomeans 1 1.2, power(0.8) alpha(0.01(0.01)0.5) sd(1) ///
	graph(y(N) x(alpha) scheme(idinsight) name(twosided,replace)) /// 
	table(_all) saving("$WD/twosided",replace) 
	
power twomeans 1 1.2, power(0.8) alpha(0.01(0.01)0.5) sd(1) onesided ///
	graph(y(N) x(alpha) scheme(idinsight) name(onesided,replace) plotop(mc(red) lcolor(red) lp(dot)) ) /// 
	table(_all) saving("$WD/onesided",replace) 
	
graph combine twosided onesided, ycommon xcommon 
graph export "$WD/TwoVSOnesidedTest.emf", replace as(emf)
