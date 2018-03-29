TITLE regular.mod   squid sodium, potassium, M-type potassium, calcium, and leak channels
 
COMMENT
 This is the regular spike without adaptation treatment based on Hodgkin-Huxley treatment
 for the set of sodium, potassium, M-type potassium, and leakage channels.
  ("Impulse encoding mechanisms of ganglion cells in the tiger salamander retina." 
  Fohlmeister.JF. 78(4):1935-47.(1997).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX regular
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
		USEION ca READ eca WRITE ica
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gabar, gcabar, gl, el, gna, gk, ga, gca
        GLOBAL minf, hinf, ninf, mtau, htau, ntau, atau, ainf, btau, binf, ctau, cinf
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .05 (S/cm2)	<0,1e9>
        gkbar = .012 (S/cm2)	<0,1e9>
		gabar = .036 (S/cm2)	<0,1e9>
		gcabar = .0022 (S/cm2)	<0,1e9>
        gl = .00005 (S/cm2)	<0,1e9>
        el = -65 (mV)
		ena = 35 (mV)
		ek = -75 (mV)
		eca = 120 (mV)
}
 
STATE {
        m h n a b c
}
 
ASSIGNED {
        v (mV)
		gna (S/cm2)
		gk (S/cm2)
		ga (S/cm2)
		gca (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
		ica (mA/cm2)
        il (mA/cm2)
        minf hinf ninf ainf binf cinf
		mtau (ms) htau (ms) ntau (ms) atau (ms) btau (ms) ctau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
		ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
		ga = gabar*a*a*a*b
		gca = gcabar*c*c*c
		ik = gk*(v - ek) + ga*(v - ek)      
        il = gl*(v - el)
		ica = gca*(v - eca)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	a = ainf
	b = binf
	c = cinf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
		a' = (ainf-a)/atau
		b' = (binf-b)/btau
		c' = (cinf-c)/ctau
}
 

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum

UNITSOFF
        alpha = .6 * vtrap(-(v+30),10)
        beta =  20 * exp(-(v+55)/18)
        sum = alpha + beta
		mtau = 1/sum
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = .4 * exp(-(v+50)/20)
        beta = 6 / (exp(-(v+20)/10) + 1)
        sum = alpha + beta
		htau = 1/sum
        hinf = alpha/sum
                :"n" potassium activation system
        alpha = .02 * vtrap(-(v+40),10) 
        beta = .4 * exp(-(v+50)/80)
		sum = alpha + beta
        ntau = 1/sum
        ninf = alpha/sum
		        :"a" A-type potassium activation system
        alpha = .006 * vtrap(-(v+90),10) 
        beta = .1 * exp(-(v+30)/10)
		sum = alpha + beta
        atau = 1/sum
        ainf = alpha/sum
		        :"b" A-type potassium activation system
        alpha = .04 * exp(-(v+70)/20)
        beta = .6 / (exp(-(v+40)/10) + 1)
		sum = alpha + beta
        btau = 1/sum
        binf = alpha/sum
		        :"c" calcium activation system
        alpha = .3 * vtrap(-(v+13),10) 
        beta = 10 * exp(-(v+38)/18)
		sum = alpha + beta
        ctau = 1/sum
        cinf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
