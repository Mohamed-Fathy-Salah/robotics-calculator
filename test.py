"""
driver file for testing
"""
from factory import *

print(">> DH parameters <<")
print("a\talpha\td\ttheta")
print()
print("0\t0\t1\tt1")
print("0\t-pi/2\td2\t0")
print("0\t0\td3\t0")

sys = Factory(0)
print("-------------- forward kinamatics --------------")
sys.fk()
input("press enter to continue")

print("-------------- inverse kinamatics --------------")
sys.ik()
input("press enter to continue")

print("-------------- jacobian --------------")
sys.jacob()
input("press enter to continue")

print("-------------- T_new --------------")
sys.j()
input("press enter to continue")

print("-------------- inverse jacobian --------------")
sys.inv_jacob()
input("press enter to continue")

print("-------------- trajectory --------------")
sys.t()
