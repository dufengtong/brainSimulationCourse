TITLE burste.mod   squid sodium, potassium, and leak channels
 
COMMENT
 This is the bursting excitatory treatment based on Hodgkin-Huxley treatment
 for the set of sodium, potassium, M-type potassium, and leakage channels.
  ("Linearization of F-I curves by adaptation" B.Ermentrout.1;10(7):1721-9.(1998).)
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
        SUFFIX burste
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gnapbar, gabar, gmbar, gnap, ga, gl, el, gna, gk, gm
        GLOBAL minf, hinf, ninf, htau, ntau, ztau, zinf, ainf, btau, binf, pinf
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .035 (S/cm2)	<0,1e9>
        gkbar = .006 (S/cm2)	<0,1e9>
		gnapbar = .0002 (S/cm2)	<0,1e9>
		gabar = .0014 (S/cm2)	<0,1e9>
        gl = .00005 (S/cm2)	<0,1e9>
		gmbar = 0.001 (S/cm2)	<0,1e9>
        el = -70 (mV)
		ena = 55 (mV)
		ek = -90 (mV)
		ztau = 75 (ms)
		btau = 15 (ms)
}
 
STATE {
        h n z b
}
 
ASSIGNED {
        v (mV)
		gna (S/cm2)
		gk (S/cm2)
		gm (S/cm2)
		gnap (S/cm2)
		ga (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf zinf ainf binf pinf
		htau (ms) ntau (ms) 
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*minf*minf*minf*h
        gk = gkbar*n*n*n*n
		gm = gmbar*z
		gnap = gnapbar*pinf
		ga = gabar*ainf*ainf*ainf*b
		ina = gna*(v - ena) + gnap*(v - ena)
		ik = gk*(v - ek) + gm*(v - ek) + ga*(v - ek)
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	h = hinf
	n = ninf
	z = zinf
	b = binf
}

? states
DERIVATIVE states {  
        rates(v)
        b' =  (binf-b)/btau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
		z' = (zinf-z)/ztau
}
 

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum

UNITSOFF
		minf = 1 / (exp(-(v+30)/9.5) + 1)
                :"h" sodium inactivation system
		htau = 0.1 + .75 / (exp((v+40.5)/6) + 1)
        hinf = 1 / (exp((v+45)/7) + 1)
                :"n" potassium activation system
        ntau = 0.1 + .5 / (exp((v+27)/15) + 1)
        ninf = 1 / (exp(-(v+35)/10) + 1)
		pinf = 1 / (exp(-(v+45)/3) + 1)
		ainf = 1 / (exp(-(v+50)/20) + 1)
		binf = 1 / (exp((v+80)/6) + 1)
		        :"z" M-type potassium activation system
        zinf = 1 / (exp(-(v+39)/5) + 1)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON