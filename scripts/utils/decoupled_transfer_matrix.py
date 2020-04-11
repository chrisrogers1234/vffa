import math
import cmath
import numpy
import copy

#citation: George Parzen 1995 "Linear Parameters and the Decoupling Matrix for
#          Linearly Coupled Motion in 6 Dimensional Phase Space",
#          arxiv:acc-phys/9510006, PAC 1997, DOI: 10.1109/PAC.1997.750716
class DecoupledTransferMatrix(object):
    """
    Calculates transformation for decoupling a 2N dimensional transfer matrix
    into N 2x2 transfer matrices.

    Once the decoupling transformation has been calculated, this can be used
    to deduce periodic solutions to the general beam ellipse problem, including
    periodic beta functions and generalised emittances.

    Member data:
    - dim is the dimension of the real space of the problem (half the size of
      the transfer matrix)
    - m is the coupled transfer matrix as defined by user, size 2*dim x 2*dim
    - det_m is the determinant of m
    - m_evalue and m_evector are the eigenvalues and corresponding eigenvectors
      of m (as calculated by numpy linalg package)
    - r is the transformation that transforms from the decoupled system to the
      coupled system
    - r_inv is the inverse of r (i.e. r^-1) that transforms from the coupled
      system to the decoupled system
    - t is the transfer matrix in the decoupled system, such that t = r^-1 m r
    - t_evalue, t_evector is the eigenvalue and corresponding eigenvectors of t
    - v_t is the periodic transfer matrix in t, such that v_t = t^T v_t t
    - phase is a list giving the phase advance in the eigensystem
    Note that from the definition of t and r, it follows that for some phase
    space vector in the coupled space u_2 = m u_1
    """
    def __init__(self, transfer_matrix):
        """
        Calculate the decoupling transformation
        - transfer_matrix: a 2Nx2N list of lists whose elements make up the
          transfer matrix. transfer_matrix should be symplectic, i.e.
          M S^T M^T S = I where S is the usual symplectic matrix.
        A ValueError is raised if det_m deviates from 1 by more than
        DecoupledTransferMatrix.det_tolerance or M is not 2N*2N. A
        FloatingPointError is raised if no periodic beta function could be
        found (e.g. lattice is not stable). 
        """
        self.m = numpy.array(transfer_matrix)
        if self.m.shape[0] != self.m.shape[1]:
            raise ValueError("Transfer matrix was not square")
        if (self.m.shape[0]/2)*2 != self.m.shape[0]:
            raise ValueError("Transfer matrix had odd size (should be even)")
        self.det_m = numpy.linalg.det(self.m)
        if abs(self.det_m - 1) > self.det_tolerance:
            raise ValueError("Transfer matrix determinant deviated from 1 by "+\
              str(self.det_m - 1)+" - DecoupledTransferMatrix.det_tolerance was "+\
              str(self.det_tolerance) )
        self.dim = int(self.m.shape[0]/2)
        self.m_evalue, self.m_evector = numpy.linalg.eig(self.m)
        self.t_evector = None
        self.t_evalue = None
        self.t = None
        self.r = None
        self.r_inv = None
        self.v_t = None
        self.phase = [None]*self.dim
        self._get_decoupling()

    def _get_decoupling(self):
        """
        Calculate the transformation that decouples the transfer matrix and some
        other associated set-up.
        """
        # par_t_evector is the t_evector from parzen paper; does not appear to
        # be consistent with numpy t_evector; I assume par_t_evector is some
        # thing just used to calculate t; but the self.t_evector is the actual
        # t_evector.
        par_t_evector = [[0+0j for i in range(self.dim*2)] for j in range(self.dim*2)]
        par_t_evector = numpy.array(par_t_evector)
        self.t_evector = numpy.array(par_t_evector)
        self.v_t = numpy.zeros([self.dim*2, self.dim*2]) # real
        for i in range(self.dim):
            j = 2*i
            #evector = numpy.transpose(self.m_evector)[i]
            evector = numpy.transpose(self.m_evector)[j]
            phi_i = numpy.angle(evector[j])
            exp_phi_i = numpy.exp(1j*phi_i)
            try:
                ratio = evector[j+1]/evector[j]
                beta_i = 1./numpy.imag(ratio)
                alpha_i = -beta_i * numpy.real(ratio)
            except FloatingPointError:
                beta_i = -1.
                alpha_i = 0.
            phase_i = numpy.angle(self.m_evalue[j])
            par_t_evector[j,   j] = cmath.sqrt(beta_i)*exp_phi_i
            par_t_evector[j+1, j] = 1./cmath.sqrt(beta_i)*(-alpha_i+1j)*exp_phi_i
            par_t_evector[j,   j+1] = numpy.conj(par_t_evector[j, j])
            par_t_evector[j+1, j+1] = numpy.conj(par_t_evector[j+1, j])
            self.phase[i] = phase_i
            sign = 1.
            if beta_i < 0:
                sign = -1.
            self.v_t[j, j] = sign*beta_i
            self.v_t[j, j+1] = -sign*alpha_i
            self.v_t[j+1, j] = -sign*alpha_i
            self.v_t[j+1, j+1] = sign*(1+alpha_i*alpha_i)/beta_i
        self.r = numpy.dot(self.m_evector, numpy.linalg.inv(par_t_evector))
        self.r_inv = numpy.linalg.inv(self.r)
        self.t = numpy.dot(self.r_inv, numpy.dot(self.m, self.r))
        self.t_evalue = numpy.array([0+0j]*(2*self.dim))
        for i in range(0, 2*self.dim, 2):
            t_quad = self.t[i:i+2, i:i+2]
            quad_evalue, quad_evector = numpy.linalg.eig(t_quad)
            self.t_evalue[i] = quad_evalue[0]
            self.t_evalue[i+1] = quad_evalue[1]
            self.t_evector[i:i+2, i:i+2] = quad_evector


    def get_v_m(self, eigen_emittances):
        """
        Get the periodic ellipse to the transfer matrix equation
        M^T V_in M = V_out, such that V_out = V_in
        - eigen_emittances: list of length N (for a transfer matrix of size
          2Nx2N). Should be filled with floats. Each float gives the 
          eigenemittance in a particular direction.
        Returns a 2Nx2N matrix V such that M^T V M = V
        """
        if len(eigen_emittances) != self.dim:
            raise ValueError("Eigen emittances "+str(eigen_emittances)+\
                    "of wrong length (should be "+str(self.dim)+")")
        v_m = copy.deepcopy(self.v_t)
        for i in range(self.dim):
            v_m[i*2, i*2] *= eigen_emittances[i]
            v_m[i*2, i*2+1] *= eigen_emittances[i]
            v_m[i*2+1, i*2] *= eigen_emittances[i]
            v_m[i*2+1, i*2+1] *= eigen_emittances[i]
        v_m = numpy.dot(numpy.dot(self.r, v_m), numpy.transpose(self.r))
        return numpy.real(v_m)


    def get_amplitude(self, axis, coupled_phase_space_vector):
        """
        Get the 2d amplitude in a given eigen plane
        """
        decoupled = self.decoupled(coupled_phase_space_vector)
        decoupled = decoupled[2*axis:2*axis+1]
        matrix = self.v_t[2*axis:2*axis+1, 2*axis:2*axis+1]
        decoupled_t = numpy.transpose(decoupled)
        amplitude = numpy.dot(numpy.dot(decoupled_t, matrix), decoupled)
        print("AMPLITUDE", amplitude.shape)
        return amplitude


    def get_phase_advance(self, axis):
        """
        Return the phase advance [rad] in the direction given by axis.
        - axis: integer, 0 >= axis < N for a transfer matrix of size 2Nx2N. 
          Defines the axis in the eigenspace along which the optics quantity is
          calculated.
        """
        return self.phase[axis]

    def get_beta(self, axis):
        """
        Return the optical beta in the direction given by axis.
        - axis: integer, 0 >= axis < N for a transfer matrix of size 2Nx2N. 
          Defines the axis in the eigenspace along which the optics quantity is
          calculated.
        """
        return self.v_t[2*axis, 2*axis]

    def get_alpha(self, axis):
        """
        Return the optical alpha in the direction given by axis.
        - axis: integer, 0 >= axis < N for a transfer matrix of size 2Nx2N. 
          Defines the axis in the eigenspace along which the optics quantity is
          calculated.
        """
        return -self.v_t[2*axis+1, 2*axis]

    def get_gamma(self, axis):
        """
        Return the optical gamma in the direction given by axis.
        - axis: integer, 0 >= axis < N for a transfer matrix of size 2Nx2N. 
          Defines the axis in the eigenspace along which the optics quantity is
          calculated.
        """
        return self.v_t[2*axis+1, 2*axis+1]

    def decoupled(self, coupled_phase_space_vector):
        """
        Transform the vector from the coupled coordinate system of M to the
        decoupled coordinate system of T
        - coupled_phase_space_vector: numpy.array() of length 2*self.dim
        Returns a numpy.array() u_m of length 2*self.dim, such that
          coupled(T decoupled(u_m) ) = M u_m
        and 
          decoupled(M coupled(u_t) ) = T u_t
        """
        coupled_phase_space_vector = numpy.array(coupled_phase_space_vector)
        decoupled_psv = numpy.dot(self.r_inv, coupled_phase_space_vector)
        return numpy.real(decoupled_psv)

    def coupled(self, decoupled_phase_space_vector):
        """
        Transform the vector u_t from the decoupled coordinate system of T to 
        the coupled coordinate system of M
        - decoupled_phase_space_vector: numpy.array() of length 2*self.dim
        Returns a numpy.array() u_m of length 2*self.dim, such that
          coupled(t decoupled(u_m) ) = m u_m
        and 
          decoupled(m coupled(u_t) ) = t u_t
        """
        decoupled_phase_space_vector = numpy.array(decoupled_phase_space_vector)
        coupled_psv = numpy.dot(self.r, decoupled_phase_space_vector)
        return numpy.real(coupled_psv)

    def print_tests(self):
        print("M")
        print(self.m)
        print("M evector")
        print(self.m_evector)
        print("M evalue")
        print(self.m_evalue)
        for i in range(self.dim):
            j = 2*i
            cos_mu = (self.t[j, j]+self.t[j+1, j+1])/2.
            print("Cos", i, cos_mu, math.cos(self.phase[i]))
            j_matrix = numpy.array([
                [self.t[j, j]-cos_mu, self.t[j, j+1]],
                [self.t[j+1, j], self.t[j+1, j+1]-cos_mu]
            ])
            det_j = cmath.sqrt(numpy.linalg.det(j_matrix))
            j_matrix /= det_j
            print("Sin", i, det_j, math.sin(self.phase[i]))
            print("Beta", i, j_matrix[0, 1], self.v_t[j, j])
            print("Alpha", i, j_matrix[0, 0], -j_matrix[1, 1], -self.v_t[j, j+1])
            print("Gamma", i, -j_matrix[1, 0], self.v_t[j+1, j+1])
        print("T evector")
        print(self.t_evector)
        print("R")
        print(self.r)
        print("T")
        print(self.t)
        print()
        print("V_T", numpy.linalg.det(self.v_t))
        print(self.v_t)
        v_t_transported = numpy.dot(self.t, numpy.dot(self.v_t, numpy.transpose(self.t)))
        print("T^T V_T T", numpy.linalg.det(v_t_transported))
        print(v_t_transported)
        print("V_M")
        v_m = self.get_v_m(list(range(2, 2+self.dim)))
        print(v_m)
        print("M^T V_M M")
        v_m_transported = numpy.dot(self.m, numpy.dot(v_m, numpy.transpose(self.m)))
        print(v_m_transported)

    @classmethod
    def symplecticity(self, matrix):
        """
        Returns matrix giving the degree of symplecticity.
        - matrix: 2N x 2N matrix
        Returns J = M S^T M^T S, where S is the symplectic matrix. For a 
        perfectly symplectic matrix this is 1.
        """
        if len(matrix.shape) != 2:
            raise ValueError("Matrix should be a matrix")
        elif matrix.shape[0] != matrix.shape[1]:
            raise ValueError("Should be a square matrix")
        elif (matrix.shape[0]/2)*2 != matrix.shape[0]:
            raise ValueError("Should be a 2N x 2N matrix")
        symp = numpy.zeros(matrix.shape)
        for i in range(1, matrix.shape[0]):
            symp[i, i-1] = -1.
            symp[i-1, i] = +1.
        matrix_T = numpy.transpose(matrix)
        # use S^T = -S
        J = numpy.dot(matrix, numpy.dot(-symp, numpy.dot(matrix_T, symp)))
        return J


    det_tolerance = 1e-6

def _random_rotated(dim):
    unit = numpy.zeros([dim, dim])
    for i in range(dim):
        unit[i, i] = 1.
    test_matrix = numpy.array(unit)
    for i, angle in enumerate(numpy.random.uniform(0., 360., dim-1)):
        rot_matrix = numpy.array(unit)
        rot_matrix[0, 0] = math.cos(angle)
        rot_matrix[i+1, i+1] = math.cos(angle)
        rot_matrix[0, i+1] = math.sin(angle)
        rot_matrix[i+1, 0] = -math.sin(angle)
        test_matrix = numpy.dot(rot_matrix, test_matrix)
    for i, angle in enumerate(numpy.random.uniform(0., 360., dim-1)):
        rot_matrix = numpy.array(unit)
        rot_matrix[0, 0] = math.cos(angle)
        rot_matrix[i+1, i+1] = math.cos(angle)
        rot_matrix[0, i+1] = math.sin(angle)
        rot_matrix[i+1, 0] = -math.sin(angle)
        test_matrix = numpy.dot(numpy.transpose(rot_matrix), test_matrix)
    return test_matrix

def test_decoupled_transport(matrix):
    """
       Check that
          coupled(T decoupled(u_m) ) = M u_m
        and 
          decoupled(M coupled(u_t) ) = T u_t
    """
    tm = DecoupledTransferMatrix(matrix)
    u_m = numpy.random.random(tm.dim*2)
    u_t = tm.decoupled(u_m)
    transported_0a = numpy.dot(tm.m, u_m)
    transported_1a = tm.coupled(numpy.dot(tm.t, tm.decoupled(u_m)))
    print("Should be equal:")
    print("    Conventional M", transported_0a)
    print("    R T R^-1 u_m  ", transported_1a)

    transported_0b = numpy.real(numpy.dot(tm.t, u_t))
    transported_1b = tm.decoupled(numpy.dot(tm.m, tm.coupled(u_t)))
    print("Should be equal:")
    print("    Conventional T", transported_0b)
    print("    R^-1 M R u_t  ", transported_1b)
    

def test_get_closed_ellipse():
    numpy.random.seed(1)
    print("============== 2D =============")
    test_matrix = [ # block-diagonal, r should be identity
    [0.75**0.5, 0.5],
    [-0.5, 0.75**0.5],
    ]
    DecoupledTransferMatrix(test_matrix).print_tests()
    test_decoupled_transport(test_matrix)
    print("\n\n============== Rotated 4D "+"="*100)
    test_matrix = [
        [0.75**0.5, 0.5,  0.0, 1.0],
        [-0.5, 0.75**0.5, -1.0, 0.0],
        [0.0, -1.0,  0.75**0.5, 0.5,],
        [1.0, 0.0,  -0.5, 0.75**0.5]
    ]
    DecoupledTransferMatrix(test_matrix).print_tests()
    test_decoupled_transport(test_matrix)
    print("\n\n============== Rotated 6D "+"="*100)
    test_matrix = _random_rotated(6)
    DecoupledTransferMatrix(test_matrix).print_tests()
    test_decoupled_transport(test_matrix)

if __name__ == "__main__":
    test_get_closed_ellipse()

