#include russian language support
import sys
reload(sys)
sys.setdefaultencoding('utf-8')
#include required libraries
from numpy import *
from pylab import *
# setup parameters and state variables
T       = 60                  
dt      = 1               	  
time    = arange(0, T+dt, dt) 
t_rest  = 0                                   
# iterate over each time step
for i, t in enumerate(time): 
  if t > t_rest:
	Vm[i] = Vm[i-1] + (-Vm[i-1] + I*Rm) / tau_m * dt
	if Vm[i] >= Vth:
	  Vm[i] += V_spike
	  t_rest = t + tau_ref