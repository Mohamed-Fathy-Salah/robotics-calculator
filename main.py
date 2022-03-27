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
    # str: di -> d1 | d2 | ...| dn ==>> d variable for prismatic joint
    def forward_kinamatics(self, joint_variables):
        # compensate in T0n with joint_variables & return x, y, z, phi, theta, psi
        tmp = self.T[-1].subs(joint_variables)
        x = tmp[0, 3]
        y = tmp[1, 3]
        z = tmp[2, 3]
        phi = atan2(tmp[1, 0], tmp[0, 0])
        psi = atan2(tmp[2, 1], tmp[2, 2])
        theta = atan2(-tmp[2,0], sqrt(1-tmp[2,0]**2))
        # print(np.array(self.T[-1].tolist()))
        return [x, y, z, phi, theta, psi]

    def inverse_kinamatics(self, x=0, y=0, z=0, phi=None, theta=None, psi=None):
        all_equations = [
                self.T[-1][0, 3] - x,
                self.T[-1][1, 3] - y,
                self.T[-1][2, 3] - z,
                ]

        if phi != None: 
            all_equations = all_equations + [
                    atan(self.T[-1][1, 0] / self.T[-1][0, 0]) - phi,
                    atan(-self.T[-1][2, 0] / sqrt(1 - self.T[-1][2, 0]**2)) - theta,
                    atan(self.T[-1][2, 1] / self.T[-1][2, 2]) - psi
                    ]

        # remove non symbolic equations 
        equations = list()
        for equation in all_equations:
            if type(equation) != float and equation !=0 and type(equation) != calculus.accumulationbounds.AccumulationBounds: 
                equations.append(equation)

        return solve(equations)


    # def jacob(self):
    #     self.__jacobian()
    #     return self.jacobian
    
    def __jacobian(self):
        if self.jacobian != None:
            return 

        self.jacobian = zeros(6,len(self.joint_type))

        O_n = self.T[-1][0:3, 3]

        if self.joint_type[0] == Joint.PRISMATIC:
            # JV_0 = z_0 , JW_0 = 0
            self.jacobian[2, 0] = 1 # z0 = [0, 0, 1]
        else:
            # JV_0 = z_0 * O_n
            self.jacobian[0:3, 0] = Matrix([0, 0, 1]).cross(O_n)

            # JW_0 = z_0
            self.jacobian[5, 0] = 1 # z0 = [0, 0, 1]

        for i in range(1, len(self.joint_type)):
            # z_i-1
            z_last = self.T[i-1][0:3, 2]
            # O_i-1
            O_last = self.T[i-1][0:3, 3]

            if self.joint_type[i] == Joint.PRISMATIC:
                # JV_i = z_i-1 , JW_i = 0
                self.jacobian[0:3, i] = z_last
            else:
                # JV_i = z_i-1 * (O_n - O_i-1)
                self.jacobian [0:3, i] = z_last.cross(O_n - O_last)

                # JW_i = z_i-1
                self.jacobian[3:6, i] = z_last

    # def inv_jacob(self):
    #     self.__inverse_jacobian()
    #     return self.inverse_jacobian
    
    def __inverse_jacobian(self):
        if self.inverse_jacobian != None:
            return
        self.__jacobian() # make sure that it is present
        self.inverse_jacobian = (self.jacobian.T * self.jacobian) ** -1 * self.jacobian.T
        
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
        self.__jacobian()
        J = np.array(self.jacobian.subs(joint_variables).tolist()).astype(np.float64)

        joint_velocity = np.array(joint_velocity).reshape(-1, 1) # assure that it is of size (n*1)

        zeta_dot = np.matmul(J , joint_velocity)
        delta = self.__delta(zeta_dot)

        t_old = self.T[-1].subs(joint_variables)
        t_new = np.add(t_old , np.dot(delta , t_old))
        return t_new

    def cubic_coeffetions(self, t0, tf, q_t0, dq_t0, q_tf, dq_tf):
        a0, a1, a2, a3 = symbols('a0 a1 a2 a3')

        solution = solve([
            Eq(q_t0, a0 + a1*t0 + a2*t0**2 + a3*t0**3),
            Eq(q_tf, a0 + a1*tf + a2*tf**2 + a3*tf**3),
            Eq(dq_t0, a1 + 2*a2*t0 + 3*a3*t0**2),
            Eq(dq_tf, a1 + 2*a2*tf + 3*a3*tf**2)
            ])
        # TODO: plot position, velocity as function of time 

    def quintic_coeffetions(self, t0, tf, q_t0, dq_t0, ddq_t0, q_tf, dq_tf, ddq_tf):
        a0, a1, a2, a3, a4, a5 = symbols('a0 a1 a2 a3 a4 a5')
        solution = solve([
            Eq(q_t0, a0 + a1*t0 + a2*t0**2 + a3*t0**3 + a4*t0**4 + a5*t0**5),
            Eq(q_tf, a0 + a1*tf + a2*tf**2 + a3*tf**3 + a4*tf**4 + a5*tf**5),
            Eq(dq_t0, a1 + 2*a2*t0 + 3*a3*t0**2 + 4*a4*t0**3 + 5*a5*t0**4),
            Eq(dq_tf, a1 + 2*a2*tf + 3*a3*tf**2 + 4*a4*tf**3 + 5*a5*tf**4),
            Eq(ddq_t0, 2*a2 + 6*a3*t0 + 12*a4*t0**2 + 20*a5*t0**3),
            Eq(ddq_tf, 2*a2 + 6*a3*tf + 12*a4*tf**2 + 20*a5*tf**3)
            ])
        # TODO: plot position, velocity, acceleration as function of time

    # end_effector_position: Matrix(x(t), y(t), z(t))
    # time: [t0, t1, ..., tn]
    def trajectory(self, end_effector_position, time):
        self.__inverse_jacobian()

        end_effector_velocity = diff(end_effector_position)

        end_effector_position_0 = np.array(end_effector_position.subs(time[0]).tolist()).astype(np.float64)
        joint_position_0 = self.inverse_kinamatics(end_effector_position_0)
        
        end_effector_velocity_0 = np.array(end_effector_velocity.subs(time[0]).tolist()).astype(np.float64)
        # TODO: multiply by zeta
        joint_velocity_0 = self.inverse_jacobian.subs(joint_position_0) * Matrix(6, 1, [])

        for i in range(1, len(time)): 
            end_effector_position_f = np.array(end_effector_position.subs(time[i]).tolist()).astype(np.float64)
            joint_position_f = self.inverse_kinamatics(end_effector_position_0)

            end_effector_velocity_f = np.array(end_effector_velocity.subs(time[0]).tolist()).astype(np.float64)
            # TODO: multiply by zeta
            joint_velocity_f = self.inverse_jacobian.subs(joint_position_f) * Matrix(6, 1, [])
            
            a0, a1, a2, a3, a4 = symbols('a0 a1 a2 a3 a4')

            for i in range(len(self.joint_type)):
                if self.joint_type[i] == Joint.PRISMATIC: 
                    pass
                else:
                    pass

            # TODO: split it for all joint variables
            tmp = solve([
                Eq(a0 + a1*time[i-1] + a2*time[i-1]**2 + a3*time[i-1]**3 + a4*time[i-1]**4, joint_position_0),
                Eq(a1 + 2*a2*time[i-1] + 3*a3*time[i-1]**2 + 4*a4*time[i-1]**3, joint_velocity_0),
                Eq(a0 + a1*time[i] + a2*time[i]**2 + a3*time[i]**3 + a4*time[i]**4, joint_position_f),
                Eq(a1 + 2*a2*time[i] + 3*a3*time[i]**2 + 4*a4*time[i]**3, joint_velocity_f)
                ])

            # TODO: plot equations 
            
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
    # system = System([DH(0, 0, 1, 0, Joint.REVOLUTE, 1),DH(0, -pi/2, 0, 0, Joint.PRISMATIC, 2), DH(0, 0, 0, 0, Joint.PRISMATIC, 3)])
    system = System([
        DH(0, 0, 1, 0, Joint.REVOLUTE, 1),
        DH(0, -pi/2, 0, 0, Joint.PRISMATIC, 2),
        DH(0, 0, 0, 0, Joint.PRISMATIC, 3)
        ])
    # system = System([
        # DH(0, -pi/2, 0, 0,Joint.REVOLUTE,1),
        # DH(0, pi/2, 2, 0,Joint.REVOLUTE,2),
        # DH(0, 0, 0, 0,Joint.PRISMATIC,3),
        # DH(0, -pi/2, 0, 0,Joint.REVOLUTE,4),
        # DH(0, pi/2, 0, 0,Joint.REVOLUTE,5),
        # DH(0, 0, 6, 0,Joint.REVOLUTE,6),
        # ])
    # system = System([
        # DH(2,0,0,.529,Joint.REVOLUTE,1),
        # DH(2,0,0,.776,Joint.REVOLUTE,2),
        # ])
    # print(system.forward_kinamatics({'t1':pi/2, 't2':pi/2, 'd3':3, 't4':pi/2, 't5':pi/2, 't6':pi/2})) 
    # print(system.inverse_kinamatics([1, -1.2, 2, 0.694738276196703, 0, 1.570796]))
    # print(system.inverse_kinamatics(x=1, y=-1.2, z=2, phi=0.694738276196703, theta=0, psi=1.570796))
    # print(system.inverse_kinamatics(x=2.5, y=1)) #(x=2.49999998504874, y=0.999999969193044, z=0, phi=-0.451835740960530, theta=0, psi=0))
    # print(system.forward_kinamatics({'t1': 1.21284847806357, 't2': -1.6646842190241}))
    # print(system.jacob())
    # print(system.inv_jacob())
    # print(system.move({'t1':pi/2, 'd2':.5, 'd3':.5}, [pi/36, .01, .01]))
    # print(system.cubic_coeffetions(0, 1, 10, 0, -20, 0))
    print(system.quintic_coeffetions(0, 2, 0, 0, 0, 40, 0, 0))

