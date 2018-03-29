TITLE regulara.mod   squid sodium, potassium, M-type potassium, and leak channels
 
COMMENT
 This is the regular spike adaptation treatment based on Hodgkin-Huxley treatment
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
        SUFFIX regulara
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gmbar, gl, el, gna, gk, gm
        GLOBAL minf, hinf, ninf, mtau, htau, ntau, ztau, zinf
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .1 (S/cm2)	<0,1e9>
        gkbar = .08 (S/cm2)	<0,1e9>
        gl = .0001 (S/cm2)	<0,1e9>
		gmbar = 0.01 (S/cm2)	<0,1e9>
        el = -67 (mV)
		ena = 50 (mV)
		ek = -100 (mV)
		ztau = 100 (ms)
}
 
STATE {
        m h n z
}
 
ASSIGNED {
        v (mV)
		gna (S/cm2)
		gk (S/cm2)
		gm (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf zinf
		mtau (ms) htau (ms) ntau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
		ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
		gm = gmbar*z
		ik = gk*(v - ek) + gm*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	z = zinf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
		z' = (zinf-z)/ztau
}
 

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum

UNITSOFF
        alpha = .32 * vtrap(-(v+54),4)
        beta =  .28 * vtrap((v+27), 5)
        sum = alpha + beta
		mtau = 1/sum
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = .128 * exp(-(v+50)/18)
        beta = 4 / (exp(-(v+27)/5) + 1)
        sum = alpha + beta
		htau = 1/sum
        hinf = alpha/sum
                :"n" potassium activation system
        alpha = .032 * vtrap(-(v+52),5) 
        beta = .5 * exp(-(v+57)/40)
		sum = alpha + beta
        ntau = 1/sum
        ninf = alpha/sum
		        :"z" M-type potassium activation system
        zinf = 1 / (exp(-(v+20)/5) + 1)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
