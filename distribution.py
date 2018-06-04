import numpy as np

class Gaussian_2D(object):

    def __init__(self, xx, xxp, xpxp):

        assert xxp**2 < xx*xpxp

        self.xx = xx
        self.xxp = xxp
        self.xpxp = xpxp

    def generate(self, n):

        sigma_matrix = np.array([[self.xx, self.xxp],[self.xxp, self.xpxp]])
        eigenvalues, eigenmatrix = np.linalg.eigh(sigma_matrix)

        self.particle_list = np.zeros( (n,2) )

        self.particle_list[:,0] = np.random.normal(scale = np.sqrt(eigenvalues[0]), size=n)
        self.particle_list[:,1] = np.random.normal(scale = np.sqrt(eigenvalues[1]), size=n)

        for i in range(n):
            self.particle_list[i] = np.dot(eigenmatrix, self.particle_list[i])

    def calculate_twiss(self):
        self.emittance = np.sqrt(self.xx*self.xpxp)
        self.beta = self.xx/self.emittance
        self.gamma = self.xpxp/self.emittance
