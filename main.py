from enum import Enum
from sympy import *
import math

class System:
    # joints    : list of Joint objects
    # DH_matrix : list of DH objects
    # A         : list of matrices contains A1, A2, ..., An
    # T         : list of matrices contains T01, T02, ...,T0n
    def __init__(self, DH_list):
        self.A = [i.get_A_matrix() for i in DH_list]

        self.T = list()
        self.T.append(self.A[0])
        for i in range(1,len(joints)):
            self.T.append(self.T[i-1] * self.A[i])

    # TODO : add other constants also ( add it from the begining in DH)
    # joint_variables: dictionary of (str ,float)
    # str: ti -> t1 | t2 | ...| tn ==>> theta variable for revolute joint
    # str: ai -> a1 | a2 | ...| an ==>> alpha variable for prismatic joint
    def forward_kinamatics(joint_variables):
        # compensate in T0n with joint_variables & return x, y, z, phi, theta, psi
        tmp = self.T[-1].subs(joint_variables)
        x = tmp[0][-1]
        y = tmp[1][-1]
        z = tmp[2][-1]
        phi = math.atan2(tmp[1][0], tmp[0][0])* 180/ math.pi
        psi = math.atan2(tmp[2][1], tmp[2][2])* 180/ math.pi
        theta = math.atan2(-tmp[2][0] * math.sin(phi), tmp[1][0])
        return [x, y, z, phi, theta, psi]

    def inverse_kinamatics():
        pass

    def jacobian():
        pass

    def trajectory():
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
            self.d = Symbol('d'+joint_number)
            self.theta = theta
        else :
            self.theta = Symbol('t'+joint_number)
            self.d = d

    def get_A_matrix():
        return Matrix([[math.cos(thata), -math.sin(theta)*math.cos(alpha), math.sin(theta)*math.sin(alpha), a*math.cos(theta)],
                [math.sin(theta), math.cos(theta)*math.sin(alpha), -math.cos(theta)*math.sin(alpha), a*math.sin(theta)],
                [0, math.sin(alpha), math.cos(alpha), d],
                [0, 0, 0, 1]])

class Joint(Enum):
    PRISMATIC = 1
    REVOLUTE = 2

if __name__ == '__main__':
    # get num of joints, joints type and DH parameters from user

    # enter 1:FK(q), 2:IK(p), 3:jacobian(), 4:trajectory()

