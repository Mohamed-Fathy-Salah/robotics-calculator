from robot import *

sys = Factory(0)
print("fk")
sys.fk()
print("-------------------")
print("ik")
sys.ik()
print("-------------------")
print("jacobian")
sys.j()
print("-------------------")
print("trajectory")
sys.t()
print("-------------------")
