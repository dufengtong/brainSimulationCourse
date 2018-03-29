TITLE fasta.mod   squid sodium, potassium, k1 potassium, k3 potassium and leak channels
 
COMMENT
 This is the fast spiking adaptation treatment based on Hodgkin-Huxley treatment
 for the set of sodium, potassium, k1 potassium, k3 potassium , and leakage channels.
  ("Function of specific K(+) channels" (1999).)
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
        SUFFIX fasta
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gk1bar, gk3bar, gl, el, gna, gk, gk1, gk3
        GLOBAL minf, hinf, ninf, mtau, htau, ntau, ptau, pinf
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .9 (S/cm2)	<0,1e9>
        gk1bar = .0018 (S/cm2)	<0,1e9>
        gl = .0041 (S/cm2)	<0,1e9>
		gk3bar = 1.8 (S/cm2)	<0,1e9>
        el = -70 (mV)
		ena = 60 (mV)
		ek = -90 (mV)
}
 
STATE {
        m h n p
}
 
ASSIGNED {
        v (mV)
		gna (S/cm2)
		gk1 (S/cm2)
		gk3 (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf pinf
		mtau (ms) htau (ms) ntau (ms) ptau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
		ina = gna*(v - ena)
        gk1 = gk1bar*n*n*n*n
		gk3 = gk3bar*p*p
		ik = gk1*(v - ek) + gk3*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	p = pinf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
		p' = (pinf-p)/ptau
}
 

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum

UNITSOFF
        alpha = 40 * vtrap(-(v-75),13.5)
        beta =  1.2262 * exp(-(v-0)/42.248)
        sum = alpha + beta
		mtau = 1/sum
        minf = alpha/sum
                :"h" sodium inactivation system
        alpha = .0035 * exp(-(v-0)/24.186)
        beta = 0.017 * vtrap(-(v+51.25),5.2)
        sum = alpha + beta
		htau = 1/sum
        hinf = alpha/sum
                :"n" potassium activation system
        alpha = .014 * vtrap(-(v+44),2.3) 
        beta = .0043 * exp(-(v+44)/34)
		sum = alpha + beta
        ntau = 1/sum
        ninf = alpha/sum
		        :"p" potassium activation system
        alpha = 1.0 * vtrap(-(v-95),11.8) 
        beta = .025 * exp(-(v-0)/22.222)
		sum = alpha + beta
        ptau = 1/sum
        pinf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
