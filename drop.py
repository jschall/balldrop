from math import e, factorial, log, exp, pi, sqrt, acosh, cosh, atan2, degrees
from numpy import sign as sgn, linspace
import pylab

g = 9.806

def v_dot(vn,ve,vd,k,wn,we):
    vappr_n = vn - wn
    vappr_e = ve - we
    vappr_d = vd
    vappr_mag = sqrt(vappr_n**2+vappr_e**2+vappr_d**2)
    if vappr_mag == 0:
        return (0.0,0.0,1.0)
    vappr_unit_n = vappr_n/vappr_mag
    vappr_unit_e = vappr_e/vappr_mag
    vappr_unit_d = vappr_d/vappr_mag
    return (-k*vappr_unit_n*vappr_mag**2,-k*vappr_unit_e*vappr_mag**2,1-k*vappr_unit_d*vappr_mag**2)

def drop_offset_numerical5(m,b,pd,target_alt,vn,ve,vd,wn,we,dt=.001):
    tolerance = .001
    global g
    h = -pd-target_alt
    k = b*h/m

    ul = h
    ut = sqrt(h/g)
    uv = sqrt(g*h)
    ua = g**2/ul

    pd = 0
    pn = 0
    pe = 0
    vn /= uv
    ve /= uv
    vd /= uv
    wn /= uv
    we /= uv
    dt /= ut

    dt_scaled = dt

    time = 0
    count=0

    while True:
        k1n,k1e,k1d = v_dot(vn,ve,vd,k,wn,we)
        k2n,k2e,k2d = v_dot(vn + (2.0/3.0)*dt_scaled*k1n, ve+(2.0/3.0)*dt_scaled*k1e,vd+(2.0/3.0)*dt_scaled*k1d, k,wn,we)

        an = (k1n+3*k2n)/4.0
        ae = (k1e+3*k2e)/4.0
        ad = (k1d+3*k2d)/4.0

        dvn = dt_scaled * an
        dve = dt_scaled * ae
        dvd = dt_scaled * ad

        vn += 0.5*dvn
        ve += 0.5*dve
        vd += 0.5*dvd

        if pd+dt_scaled*vd >= 1:
            dt_scaled = (1.0-pd)/vd

        time += dt_scaled

        pn += dt_scaled*vn
        pe += dt_scaled*ve
        pd += dt_scaled*vd

        #print pd, vn, ve, vd

        vn += 0.5*dvn
        ve += 0.5*dve
        vd += 0.5*dvd
        count+=1
        
        if pd >= 1:
            print count
            return (pn*ul, pe*ul, time*ut)

#ball constants
#mass of a tennis ball is 58 grams
m = 0.0568
#Fd = b*v^2, let b = 0.5*Cd*p*A
b = 0.5*.531*1.25*.003525
#drop conditions
vn = 10
ve = 0.0
vd = 0.0
wn = 0.0
we = 0.0
pd = -10000

offset, _,time = drop_offset_numerical5(m,b,pd,0.0,vn,ve,vd,wn,we,.1)
#offset, _, time = drop_offset_numerical5(m,b,pd,0.0,vn,ve,vd,wn,we,.1)

print offset, time

#tana = []
#tsim = []
#pana = []
#psim = []
#hana = []
#tana = []
#h = []



#for pd in linspace(0,-100,100):
    #px_simulated,_,time_simulated =    drop_offset_simulated(m,b,pd,0.0,vn,ve,vd,wn,we)
    #px_analytical,_,time_analytical = drop_offset_analytical(m,b,pd,0.0,vn,ve,vd,wn,we)
    ##print pd, px_simulated

    #tsim.append(time_simulated)
    #pana.append(px_analytical)
    #psim.append(px_simulated)
    #tana.append(time_analytical)
    #h.append(-pd)

#pylab.plot(h,psim,'r',h,pana,'b')
#pylab.show()
