import numpy as np
from enum import Enum
from sympy import *
from sympy.plotting import PlotGrid

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
        
        self.jacob = None
        self.inv_jacob = None

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
        print("x:{} y:{} z:{} phi:{} theta:{} psi:{}".format(x, y, z, phi, theta, psi))
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

        solutions =  solve(equations, dict=True)
        print(solutions)
        return solutions 

    def jacobian(self):
        if self.jacob!= None:
            return self.jacob

        self.jacob = zeros(6,len(self.joint_type))

        O_n = self.T[-1][0:3, 3]

        if self.joint_type[0] == Joint.PRISMATIC:
            # JV_0 = z_0 , JW_0 = 0
            self.jacob[2, 0] = 1 # z0 = [0, 0, 1]
        else:
            # JV_0 = z_0 * O_n
            self.jacob[0:3, 0] = Matrix([0, 0, 1]).cross(O_n)

            # JW_0 = z_0
            self.jacob[5, 0] = 1 # z0 = [0, 0, 1]

        for i in range(1, len(self.joint_type)):
            # z_i-1
            z_last = self.T[i-1][0:3, 2]
            # O_i-1
            O_last = self.T[i-1][0:3, 3]

            if self.joint_type[i] == Joint.PRISMATIC:
                # JV_i = z_i-1 , JW_i = 0
                self.jacob[0:3, i] = z_last
            else:
                # JV_i = z_i-1 * (O_n - O_i-1)
                self.jacob[0:3, i] = z_last.cross(O_n - O_last)

                # JW_i = z_i-1
                self.jacob[3:6, i] = z_last
        print(np.array(self.jacob.tolist()))
        return self.jacob

    def inverse_jacobian(self):
        if self.inv_jacob != None:
            return self.inv_jacob
        self.jacobian() # make sure that it is present
        self.inv_jacob = (self.jacob.T * self.jacob) ** -1 * self.jacob.T
        print(np.array(self.inv_jacob))
        return self.inv_jacob
        
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
        self.jacobian()
        J = np.array(self.jacob.subs(joint_variables).tolist()).astype(np.float64)
        print("jacobian: \n{}".format(J))

        joint_velocity = np.array(joint_velocity).reshape(-1, 1) # assure that it is of size (n*1)

        zeta_dot = np.matmul(J , joint_velocity)
        delta = self.__delta(zeta_dot)
        print("delta: \n{}".format(delta))

        t_old = self.T[-1].subs(joint_variables)
        print("t_old: \n{}".format(np.array(t_old.tolist())))

        t_new = np.add(t_old , np.dot(delta , t_old))
        print("t_new: \n{}".format(t_new))
        return t_new

    def cubic_coeffetions(self, t0, tf, q_t0, dq_t0, q_tf, dq_tf, joint_number=1):
        a0, a1, a2, a3, t = symbols('a0 a1 a2 a3 t')
        solution = solve([
            Eq(q_t0, a0 + a1*t0 + a2*t0**2 + a3*t0**3),
            Eq(q_tf, a0 + a1*tf + a2*tf**2 + a3*tf**3),
            Eq(dq_t0, a1 + 2*a2*t0 + 3*a3*t0**2),
            Eq(dq_tf, a1 + 2*a2*tf + 3*a3*tf**2)
            ], dict=True)[0]
        position = solution[a0] + solution[a1]*t + solution[a2]*t**2 + solution[a3]*t**3
        velocity = solution[a1] + 2*solution[a2]*t + 3*solution[a3]*t**2
        acceleration = 2*solution[a2] + 6*solution[a3]*t
        init_printing()
        print("\n-- joint{}".format(joint_number))
        print("position = {}".format(position))
        print("velocity = {}".format(velocity))
        print("acceleration = {}".format(acceleration))
        return [
                plot(position, show=False, title="joint{} position".format(joint_number)),
                plot(velocity, show=False, title="joint{} velocity".format(joint_number)),
                plot(acceleration, show=False, title="joint{} acceleration".format(joint_number))]
                

    def quintic_coeffetions(self, t0, tf, q_t0, dq_t0, ddq_t0, q_tf, dq_tf, ddq_tf):
        a0, a1, a2, a3, a4, a5 = symbols('a0 a1 a2 a3 a4 a5')
        solution = solve([
            Eq(q_t0, a0 + a1*t0 + a2*t0**2 + a3*t0**3 + a4*t0**4 + a5*t0**5),
            Eq(q_tf, a0 + a1*tf + a2*tf**2 + a3*tf**3 + a4*tf**4 + a5*tf**5),
            Eq(dq_t0, a1 + 2*a2*t0 + 3*a3*t0**2 + 4*a4*t0**3 + 5*a5*t0**4),
            Eq(dq_tf, a1 + 2*a2*tf + 3*a3*tf**2 + 4*a4*tf**3 + 5*a5*tf**4),
            Eq(ddq_t0, 2*a2 + 6*a3*t0 + 12*a4*t0**2 + 20*a5*t0**3),
            Eq(ddq_tf, 2*a2 + 6*a3*tf + 12*a4*tf**2 + 20*a5*tf**3)
            ], dict=True)[0]
        return solution

    # end_effector_position: Matrix(x(t), y(t), z(t))
    # time: [t0, t1, ..., tn]
    def trajectory(self, end_effector_position, time):
        self.inverse_jacobian()

        end_effector_velocity = diff(end_effector_position)

        # x(t0), y(t0), z(t0)
        end_effector_position_0 = np.array(end_effector_position.subs({'t': time[0]}).tolist()).astype(np.float64).reshape(-1,) 
        # q(t0) : list(dictionary)
        joint_position_0 = self.inverse_kinamatics(end_effector_position_0[0], end_effector_position_0[1], end_effector_position_0[2])[-1]
        
        # dx(t0), dy(t0), dz(t0)
        end_effector_velocity_0 = np.array(end_effector_velocity.subs({'t': time[0]}).tolist()).astype(np.float64).reshape(1,-1)
        # dq(t0) : n*1 array
        joint_velocity_0 = self.inv_jacob.subs(joint_position_0) * Matrix(6, 1, np.append(end_effector_velocity_0, [0, 0, 0]))

        for t in range(1, len(time)): 
            # x(tf), y(tf), z(tf)
            end_effector_position_f = np.array(end_effector_position.subs({'t': time[t]}).tolist()).astype(np.float64).reshape(-1,)
            # q(tf) : list(dictionary)
            joint_position_f = self.inverse_kinamatics(end_effector_position_f[0], end_effector_position_f[1], end_effector_position_f[2])[-1]

            # dx(tf), dy(tf), dz(tf)
            end_effector_velocity_f = np.array(end_effector_velocity.subs({'t': time[t]}).tolist()).astype(np.float64).reshape(1,-1)
            # dq(tf) : n*1 array
            joint_velocity_f = self.inv_jacob.subs(joint_position_f) * Matrix(6, 1, np.append(end_effector_velocity_f, [0, 0, 0]))
            
            # a0, a1, a2, a3, a4 = symbols('a0 a1 a2 a3 a4')
            print('\n from time t_{}={} to time t_{}={}'.format(t-1, time[t-1], t, time[t]))
            # joint_equations = None
            joint_equations = list()
            for j in range(len(self.joint_type)):
                q = None
                if self.joint_type[j] == Joint.PRISMATIC: 
                    q = Symbol('d{}'.format(j+1))
                else:
                    q = Symbol('t{}'.format(j+1))
                
                joint_equations = joint_equations + self.cubic_coeffetions(time[t-1], time[t], joint_position_0[q], joint_velocity_0[j], joint_position_f[q], joint_velocity_f[j], j+1)

                # print('\n\t>>> joint_{}<<<\n'.format(j+1))
                # print(time[t-1], time[t], joint_position_0[q] * mul, joint_velocity_0[j], joint_position_f[q] * mul, joint_velocity_f[j])
                # print('\t{}'.format(solution))
                
                # TODO: print the table

            plt = PlotGrid(1, len(joint_equations), *joint_equations)

            joint_position_0 = joint_position_f  
            joint_velocity_0 = joint_velocity_f

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
    # print(system.inverse_kinamatics(x=2.5, y=1)) #(x=2.49999998504874, y=0.999999969193044, z=0, phi=-0.451835740960530, theta=0, psi=0))
    # print(system.jacob())
    # print(system.inv_jacob())
    pass
