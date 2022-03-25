import numpy as np
from enum import Enum
from sympy import *

class System:
    # DH_matrix  : list of DH objects
    # A          : list of matrices contains A1, A2, ..., An
    # T          : list of matrices contains T01, T02, ..., T0n
    # joint_type : list of Joint object
    def __init__(self, DH_list):
        self.joint_type = [i.joint_type for i in DH_list]

        self.A = [i.get_A_matrix() for i in DH_list]

        self.T = list()
        self.T.append(self.A[0])
        for i in range(1, len(DH_list)):
            self.T.append(self.T[i-1] * self.A[i])
        
        self.jacobian = None
        self.inverse_jacobian = None

    # joint_variables: dictionary of (str, float)
    # str: ti -> t1 | t2 | ...| tn ==>> theta variable for revolute joint
    # str: ai -> a1 | a2 | ...| an ==>> alpha variable for prismatic joint
    def forward_kinamatics(self, joint_variables):
        # compensate in T0n with joint_variables & return x, y, z, phi, theta, psi
        tmp = self.T[-1].subs(joint_variables)
        x = tmp[0, 3]
        y = tmp[1, 3]
        z = tmp[2, 3]
        phi = atan2(tmp[1, 0], tmp[0, 0])* 180/ pi
        psi = atan2(tmp[2, 1], tmp[2, 2])* 180/ pi
        theta = atan2(-tmp[2, 0] * sin(phi), tmp[1, 0])
        return [x, y, z, phi, theta, psi]

    # TODO: fix it does not get all acceptable outputs maybe use another function
    # end_effector: [x, y, z, phi, theta, psi]
    def inverse_kinamatics(self, end_effector): 
        if len(end_effector) == 3:
            return solve([
                Eq(end_effector[0], self.T[-1][0, 3]),
                Eq(end_effector[1], self.T[-1][1, 3]),
                Eq(end_effector[2], self.T[-1][2, 3])
                ])

        return solve([
            Eq(end_effector[0], self.T[-1][0, 3]),
            Eq(end_effector[1], self.T[-1][1, 3]),
            Eq(end_effector[2], self.T[-1][2, 3]),
            Eq(tan(end_effector[3]), self.T[-1][1, 0]/self.T[-1][0, 0]),
            Eq(sin(end_effector[4]), -self.T[-1][3, 0]),
            Eq(tan(end_effector[5]), self.T[-1][2, 1]/self.T[-1][2, 2])
            ])

    # def jacob(self, joint_variables):
    #     return self.__jacobian(joint_variables)

    # TODO: make jacobian only once during the life time of the program
    def __jacobian(self, joint_variables):
        if self.jacobian != None:
            return 

        self.jacobian = np.zeros((6, len(self.joint_type)))
        T_n = self.T[-1].subs(joint_variables)
        O_n = [T_n[0, 3], T_n[1, 3], T_n[2, 3]]

        if self.joint_type[0] == Joint.PRISMATIC:
            # JV_0 = z_0 , JW_0 = 0
            self.jacobian[2][0] = 1 # z0 = [0, 0, 1]
        else:
            # JV_0 = z_0 * O_n
            self.jacobian [0:3, 0] = np.cross(np.array([0, 0, 1]), np.array(O_n))

            # JW_0 = z_0
            self.jacobian[5][0] = 1 # z0 = [0, 0, 1]

        for i in range(1, len(self.joint_type)):
            # z_i-1
            tmp = self.T[i-1].subs(joint_variables)
            z_last = [tmp[0, 2], tmp[1, 2], tmp[2, 2]]
            O_last = [tmp[0, 3], tmp[1, 3], tmp[2, 3]]

            if self.joint_type[i] == Joint.PRISMATIC:
                # JV_i = z_i-1 , JW_i = 0
                self.jacobian[0:3, i] = z_last
            else:
                # JV_i = z_i-1 * (O_n - O_i-1)
                self.jacobian [0:3, i] = np.cross(np.array(z_last), np.array(O_n - O_last))

                # JW_i = z_i-1
                self.jacobian[3:6, i] = z_last
    
    # zeta_dot: list of shape(6*1)
    def __delta(self, zeta_dot):
        return np.array([
            [0, float(-zeta_dot[5]), float(zeta_dot[4]), float(zeta_dot[0])],
            [float(zeta_dot[5]), 0, float(-zeta_dot[3]), float(zeta_dot[1])],
            [float(-zeta_dot[4]), float(zeta_dot[3]), 0, float(zeta_dot[2])],
            [0, 0, 0, 0] ])

    # joint_variables: dictionary of (str, float)
    # joint_velocity : list of shape (n*1)
    def move(self, joint_variables, joint_velocity):
        self.__jacobian(joint_variables)

        joint_velocity = np.array(joint_velocity).reshape(-1, 1) # assure that it is of size (n*1)

        zeta_dot = np.matmul(jacobian , joint_velocity)
        delta = self.__delta(zeta_dot)

        t_old = self.T[-1].subs(joint_variables)
        t_new = np.add(t_old , np.dot(delta , t_old))
        return t_new

    def trajectory(self):
        pass

class DH:
    # a, alpha, d, theta: float
    # joint_type: Joint
    # joint_number: int
    def __init__(self, a, alpha, d, theta, joint_type, joint_number):
        self.a = a
        self.alpha = alpha
        self.joint_type = joint_type
        # self.joint_number = joint_number
        if joint_type == Joint.PRISMATIC :
            self.d = Symbol('d{}'.format(joint_number)) + d
            self.theta = theta
        else :
            self.d = d
            self.theta = Symbol('t{}'.format(joint_number))

    def get_A_matrix(self):
        return Matrix([[cos(self.theta), -sin(self.theta)*cos(self.alpha), sin(self.theta)*sin(self.alpha), self.a*cos(self.theta)],
                [sin(self.theta), cos(self.theta)*cos(self.alpha), -cos(self.theta)*sin(self.alpha), self.a*sin(self.theta)],
                [0.0, sin(self.alpha), cos(self.alpha), self.d],
                [0.0, 0.0, 0.0, 1.0]])

class Joint(Enum):
    PRISMATIC = 1
    REVOLUTE = 2

if __name__ == '__main__':
    # get num of joints, joints type and DH parameters from user
    # num_of_joints = int(input("enter number of joints: "))

    # DH_list = list()
    # for i in range(num_of_joints):
        # joint_type = int(input("enter joint {} type [1:prismatic|2:revolute]".format(i+1)))
        # print(joint_type)
        # a = float(input("a{} :".format(i+1)))
        # alpha = float(input("alpha{} :".format(i+1)))
        # d = float(input("d{} :".format(i+1)))
        # if(joint_type == Joint.PRISMATIC.value):
            # theta = float(input("theta{} :".format(i+1)))
        # DH_list.append(DH(a, alpha, d, theta, joint_type, i))
    
    # enter 1:FK(q), 2:IK(p), 3:jacobian(), 4:trajectory()
    # system = System(DH_list)
    system = System([DH(0, 0, 1, 0, Joint.REVOLUTE, 1),DH(0, -pi/2, 0, 0, Joint.PRISMATIC, 2), DH(0, 0, 0, 0, Joint.PRISMATIC, 3)])
    # print(system.forward_kinamatics({'t1':pi/2, 'd2':.1, 'd3':.1})) 
    print(system.inverse_kinamatics([1, -1.2, 2]))
    # print(system.jacob({'t1':pi/2, 'd2':.5, 'd3':.5}))
    # print(system.move({'t1':pi/2, 'd2':.5, 'd3':.5}, [pi/36, .01, .01]))
    
