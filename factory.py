"""
this class outputs robot samples ( for testing )
use :
    system = Factory(0)
    system.fk()
    system.ik()
    system.jacob()
    system.inv_jacob()
    system.j()
    system.t()
"""
from robot import *
class Factory:
    def __init__(self, i):
        self.i = i
        self.system = self.get_system()

    def get_system(self):
        return [
                System([
                    DH(0, 0, 1, 0, Joint.REVOLUTE, 1),
                    DH(0, -pi/2, 0, 0, Joint.PRISMATIC, 2),
                    DH(0, 0, 0, 0, Joint.PRISMATIC, 3)
                    ]),
                System([
                    DH(0.1, 0, .2, 0 , Joint.REVOLUTE, 1),
                    DH(0.2, pi / 2, 0, 0, Joint.PRISMATIC, 2),
                ]),
                System([
                    DH(0, -pi/2, 0, 0,Joint.REVOLUTE,1),
                    DH(0, pi/2, 2, 0,Joint.REVOLUTE,2),
                    DH(0, 0, 0, 0,Joint.PRISMATIC,3),
                    DH(0, -pi/2, 0, 0,Joint.REVOLUTE,4),
                    DH(0, pi/2, 0, 0,Joint.REVOLUTE,5),
                    DH(0, 0, 6, 0,Joint.REVOLUTE,6),
                    ]),
                System([
                    DH(2,0,0,0,Joint.REVOLUTE,1),
                    DH(2,0,0,0,Joint.REVOLUTE,2),
                    ])
                ][self.i]

    def fk(self):
        return self.system.forward_kinamatics(self.get_joint_variables())

    def ik(self):
        return self.system.inverse_kinamatics(*self.get_end_effector())
    
    def jacob(self):
        return self.system.jacobian()

    def inv_jacob(self):
        return self.system.inverse_jacobian()

    def j(self):
        return self.system.move(self.get_joint_variables(), self.get_joint_velocity())

    def t(self):
        return self.system.trajectory(self.get_end_effector_equaitons(), self.get_time())

    def get_joint_variables(self):
        return [
                {'t1':pi/2, 'd2':.5, 'd3':.5},
                {'t1':pi/2, 'd2':.3},
                {'t1':pi/2, 't2':pi/2, 'd3':3, 't4':pi/2, 't5':pi/2, 't6':pi/2},
                {'t1': 1.21284847806357, 't2': -1.6646842190241}
                ][self.i]

    def get_end_effector(self):
        #TODO: add all 3 
        return [
                [1, -1.2, 2, 0.694738276196703, 0, 1.570796],
                [0, .3, .5, pi/2, 0, pi/2],
                [1, -1.2, 2, 0, pi/2, 0],
                ][self.i]

    def get_cubic_vars(self):
        return [0, 1, 10, 0, -20, 0]

    def get_quintic_vars(self):
        return [0, 2, 0, 0, 0, 40, 0, 0]

    def get_joint_velocity(self):
        #TODO: add all 3
        return [
                [0, pi/8, pi/4],
                ][self.i]

    def get_end_effector_equaitons(self):
        t = Symbol('t')
        return [
                Matrix([200*t, 640-160*t*t, 600]),
                Matrix([2+.5*cos(t), 1+.5*sin(t)])
                ][self.i]

    def get_time(self):
        return [
                [0, .25, .5],
                [0, pi/8, pi/4]
                ][self.i]
