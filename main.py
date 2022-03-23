from enum import Enum
from sympy import *

class System:
    # DH_matrix : list of DH objects
    # A         : list of matrices contains A1, A2, ..., An
    # T         : list of matrices contains T01, T02, ...,T0n
    def __init__(self, DH_list):
        self.A = [i.get_A_matrix() for i in DH_list]
        self.T = list()
        self.T.append(self.A[0])
        for i in range(1,len(DH_list)):
            self.T.append(self.T[i-1] * self.A[i])

    # TODO : add other constants also ( add it from the begining in DH)
    # joint_variables: dictionary of (str ,float)
    # str: ti -> t1 | t2 | ...| tn ==>> theta variable for revolute joint
    # str: ai -> a1 | a2 | ...| an ==>> alpha variable for prismatic joint
    def forward_kinamatics(self,joint_variables):
        # compensate in T0n with joint_variables & return x, y, z, phi, theta, psi
        tmp = self.T[-1].subs(joint_variables)
        x = tmp[0,3]
        y = tmp[1,3]
        z = tmp[2,3]
        phi = atan2(tmp[1,0], tmp[0,0])* 180/ pi
        psi = atan2(tmp[2,1], tmp[2,2])* 180/ pi
        theta = atan2(-tmp[2,0] * sin(phi), tmp[1,0])
        return [x, y, z, phi, theta, psi]

    def inverse_kinamatics(self):
        pass

    def jacobian(self):
        pass

    def trajectory(self):
        pass

class DH:
    # a, alpha, d, theta: float
    # joint_type: Joint
    # joint_number: int
    def __init__(self, a, alpha, d, theta, joint_type, joint_number):
        self.a = a
        self.alpha = alpha
        self.joint_type = joint_type # TODO: not needed ? 
        self.joint_number = joint_number # TODO: not needed ?
        if joint_type == Joint.PRISMATIC :
            self.d = Symbol('d{}'.format(joint_number))
            self.theta = theta
        else :
            self.theta = Symbol('t{}'.format(joint_number))
            self.d = d

    def get_A_matrix(self):
        # return Matrix([[1,0,0,0],[0,self.theta,0,0],[0,0,self.theta,0],[0,0,0,self.theta]])
        return Matrix([[cos(self.theta), -sin(self.theta)*cos(self.alpha), sin(self.theta)*sin(self.alpha), self.a*cos(self.theta)],
                [sin(self.theta), cos(self.theta)*sin(self.alpha), -cos(self.theta)*sin(self.alpha), self.a*sin(self.theta)],
                [0.0, sin(self.alpha), cos(self.alpha), self.d],
                [0.0, 0.0, 0.0, 1.0]])

class Joint(Enum):
    PRISMATIC = 1
    REVOLUTE = 2

if __name__ == '__main__':
    # get num of joints, joints type and DH parameters from user
    num_of_joints = int(input("enter number of joints: "))

    DH_list = list()
    for i in range(num_of_joints):
        joint_type = int(input("enter joint {} type [1:prismatic|2:revolute]".format(i+1)))
        print(joint_type)
        a = float(input("a{} :".format(i+1)))
        alpha = float(input("alpha{} :".format(i+1)))
        if(joint_type == Joint.PRISMATIC.value):
            theta = float(input("theta{} :".format(i+1)))
            d = None
        else :
            d = float(input("d{} :".format(i+1)))
            theta = None
        DH_list.append(DH(a,alpha,d,theta,joint_type,i))
    
    # enter 1:FK(q), 2:IK(p), 3:jacobian(), 4:trajectory()
    system = System(DH_list)
    # system = System([DH(0,0,.5,0,Joint.REVOLUTE,1),DH(0,-90,0,0,Joint.PRISMATIC,2),DH(0,0,0,0,Joint.PRISMATIC,3)])
    # print(system.forward_kinamatics({'t1':pi/2,'d2':.1,'d3':.1}))
    
