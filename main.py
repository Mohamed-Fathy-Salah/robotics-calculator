from dataclasses import dataclass
from enum import Enum
from sympy import *

class System:
    # joints    : list of Joint objects
    # DH_matrix : list of DH objects
    # A         : list of matrices contains A1, A2, ..., An
    # T         : list of matrices contains T01, T02, ...,T0n
    def __init__(self, joints, DH_list):
        assert len(joints) == len(DH_list)

        self.joints = joints

        self.A = [i.get_A_matrix() for i in DH_list]

        self.T = list()
        self.T.append(self.A[0])
        for i in range(1,len(joints)):
            self.T.append(self.T[i-1] * self.A[i])

    # TODO : add other constants also ( add it from the begining in DH)
    def forward_kinamatics(joint_variables):
        # compensate in T0n with joint_variables & return x, y, z, phi, theta, psi
        pass

    def inverse_kinamatics():
        pass

    def jacobian():
        pass

    def trajectory():
        pass


@dataclass
class DH:
    a: str
    alpha: str
    d: str
    theta: str
    def get_A_matrix():
        # TODO : generate A matrix for this joint
        pass

class Joint(Enum):
    PRISMATIC = 1
    REVOLUTE = 2

