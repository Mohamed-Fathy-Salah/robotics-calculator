"""
driver file for testing
"""
from factory import Factory

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
