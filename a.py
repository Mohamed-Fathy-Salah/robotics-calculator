from sympy import *
from robot import *

system = System([
    DH(2,0,0,0,Joint.REVOLUTE,1),
    DH(2,0,0,0,Joint.REVOLUTE,2)
    ])

print(np.array(system.jacobian().subs({'t1':-20.4*pi/180,'t2':84.6*pi/180}).tolist()).astype(float))
print(np.array(system.inverse_jacobian().subs({'t1':-20.4*pi/180,'t2':84.6*pi/180}).tolist()).astype(float))
