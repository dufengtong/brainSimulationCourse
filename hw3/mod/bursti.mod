TITLE bursti.mod   squid sodium, potassium, and leak channels
 
COMMENT
 This is the bursting inhibitory treatment based on Hodgkin-Huxley treatment
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
        SUFFIX bursti
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gabar, gdbar, gd, gl, el, gna, gk
        GLOBAL minf, hinf, ninf, htau, ntau, atau, ainf, binf, btau
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .1125 (S/cm2)	<0,1e9>
        gkbar = .225 (S/cm2)	<0,1e9>
		gdbar = .0025 (S/cm2)	<0,1e9>
        gl = .00025 (S/cm2)	<0,1e9>
        el = -70 (mV)
		ena = 50 (mV)
		ek = -90 (mV)
		atau = 2 (ms)
		btau = 150 (ms)
}
 
STATE {
        h n a b
}
 
ASSIGNED {
        v (mV)
		gna (S/cm2)
		gk (S/cm2)
		gd (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf ainf binf 
		htau (ms) ntau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*minf*minf*minf*h
        gk = gkbar*n*n
		gd = gdbar*a*a*a*b
		ina = gna*(v - ena) 
		ik = gk*(v - ek) + gd*(v - ek) 
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	h = hinf
	n = ninf
	a = ainf
	b = binf
}

? states
DERIVATIVE states {  
        rates(v)
        b' =  (binf-b)/btau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
		a' = (ainf-a)/atau
}
 

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum

UNITSOFF
		minf = 1 / (exp(-(v+24)/11.5) + 1)
                :"h" sodium inactivation system
		htau = 0.5 + 14 / (exp((v+60)/12) + 1)
        hinf = 1 / (exp((v+58.3)/6.7) + 1)
                :"n" potassium activation system
        ntau = (0.087 + 11.4 / (exp((v+14.6)/8.6) + 1)) * (0.087 + 11.4 / (exp(-(v-1.3)/18.7) + 1))
        ninf = 1 / (exp(-(v+12.4)/6.8) + 1)
		ainf = 1 / (exp(-(v+50)/20) + 1)
		binf = 1 / (exp((v+70)/6) + 1)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON