import numpy as np
import matplotlib.pyplot as plt

class invariant_ellipse(object):
    
    def __init__(self, matrix):
        
        self.transfer_matrix = matrix
        
        self.mu = np.arccos(np.trace(matrix)/2.)
        self.alpha = (matrix[0,0] - np.cos(self.mu))/np.sin(self.mu)
        self.beta = matrix[0,1]/np.sin(self.mu)
        self.gamma = (1. + self.alpha**2)/self.beta 

        self.emittance = 1.
        self.calculate_moments()
        
        self.ellipse_matrix = np.array([[self.gamma, self.alpha],[self.alpha, self.beta]])     
        
        self.generate_test()
        self.generate_mismatched()
        self.generate_matched()
    
        self.generate_inv_ellipse()
        self.generate_effective_ellipse()
    
    def generate_test(self):
        self.test_particles_list = np.array([[np.sqrt(self.emittance/self.gamma), 0],
                                             [-np.sqrt(self.emittance/self.gamma), 0],
                                             [0, np.sqrt(self.emittance/self.beta)],
                                             [0, -np.sqrt(self.emittance/self.beta)]])
    
    def generate_mismatched(self):
        
        self.mismatched_particles_list = self.generate(xx=self.xpxp, xxp=-self.xxp, xpxp=self.xx, n=1000)
    
    def generate_matched(self):
        
        self.matched_particles_list = self.generate(xx=self.xx, xxp=self.xxp, xpxp=self.xpxp, n=1000)
    
    def calculate_moments(self):
        
        self.xx = self.beta*self.emittance
        self.xpxp = self.gamma*self.emittance
        self.xxp = np.sqrt(self.xx*self.xpxp - self.emittance)
        
    def generate(self, xx, xxp, xpxp, n):
        
        sigma_matrix = np.array([[xx, xxp],[xxp, xpxp]])
        eigenvalues, eigenmatrix = np.linalg.eigh(sigma_matrix)

        r1 = np.sqrt(eigenvalues[0])
        r2 = np.sqrt(eigenvalues[1])
        
        particle_list = np.zeros( (n,2) )

        i = 0    
    
        while i < n:
            x1 = np.random.uniform(-1.0,1.0)
            x2 = np.random.uniform(-1.0,1.0)
            S = (x1**2+x2**2)**0.5
            
            if S <= 1:
                particle_list[i,0] = x1*r1
                particle_list[i,1] = x2*r2
                
                i += 1
                
        for i in range(n):
            particle_list[i] = np.dot(eigenmatrix, particle_list[i])
            
        return particle_list
    
    def generate_inv_ellipse(self):

        eigen_values, eigen_vector_matrix = np.linalg.eig(self.ellipse_matrix)
        rotation_matrix = eigen_vector_matrix

        circle_list = np.array([ [np.cos(theta), np.sin(theta)] for theta in np.linspace(0., 2*np.pi, 201) ])
        ellipse_u = np.array([[u0[0]*np.sqrt(self.emittance/eigen_values[0]),
                              u0[1]*np.sqrt(self.emittance/eigen_values[1])] for u0 in circle_list])
        self.ellipse_x = np.array( [ np.dot( rotation_matrix , u)for u in ellipse_u ] )
    
    def generate_effective_ellipse(self):
        
        eigen_values, eigen_vector_matrix = np.linalg.eig(self.ellipse_matrix)
        rotation_matrix = eigen_vector_matrix

        effective_emittance = self.emittance*(np.amax(eigen_values)/np.amin(eigen_values))
        
        circle_list = np.array([ [np.cos(theta), np.sin(theta)] for theta in np.linspace(0., 2*np.pi, 201) ])
        ellipse_u = np.array([[u0[0]*np.sqrt(effective_emittance/eigen_values[0]),
                              u0[1]*np.sqrt(effective_emittance/eigen_values[1])] for u0 in circle_list])
        self.effective_ellipse_x = np.array( [ np.dot( rotation_matrix , u)for u in ellipse_u ] )        
    
    def plot_test_particle(self):
        
        plt.plot(self.ellipse_x[:,0], self.ellipse_x[:,1])
        plt.scatter(self.test_particles_list[0,0], self.test_particles_list[0,1], color='r')
        plt.scatter(self.test_particles_list[1,0], self.test_particles_list[1,1], color='k')
        plt.scatter(self.test_particles_list[2,0], self.test_particles_list[2,1], color='g')
        plt.scatter(self.test_particles_list[3,0], self.test_particles_list[3,1], color='y')
        plt.xlabel('x')
        plt.ylabel("x'")
        plt.show()        
    
    def plot_mismatched(self):
        
        plt.plot(self.ellipse_x[:,0], self.ellipse_x[:,1])
        plt.plot(self.effective_ellipse_x[:,0], self.effective_ellipse_x[:,1])
        plt.scatter(self.mismatched_particles_list[:,0], self.mismatched_particles_list[:,1], color='r', s=0.5)
        plt.xlabel('x')
        plt.ylabel("x'")
        plt.show()
    
    def plot_matched(self):
        
        plt.plot(self.ellipse_x[:,0], self.ellipse_x[:,1])
        plt.scatter(self.matched_particles_list[:,0], self.matched_particles_list[:,1], color='r', s=0.5)
        plt.xlabel('x')
        plt.ylabel("x'")
        plt.show()
        
    def plot_both(self):
        
        fig = plt.figure(figsize=(14,6))

        ax1 = fig.add_subplot(121)       
        ax1.plot(self.ellipse_x[:,0], self.ellipse_x[:,1])
        ax1.plot(self.effective_ellipse_x[:,0], self.effective_ellipse_x[:,1])
        ax1.scatter(self.matched_particles_list[:,0], self.matched_particles_list[:,1], color='r', s=0.5)
        ax1.set_xlabel('x')
        ax1.set_ylabel("x'")
        ax1.set_title('matched')

        ax2 = fig.add_subplot(122)
        ax2.plot(self.ellipse_x[:,0], self.ellipse_x[:,1])
        ax2.plot(self.effective_ellipse_x[:,0], self.effective_ellipse_x[:,1])
        ax2.scatter(self.mismatched_particles_list[:,0], self.mismatched_particles_list[:,1], color='r', s=0.5)
        ax2.set_xlabel('x')
        ax2.set_ylabel("x'")
        ax2.set_title('mismatched')
        
        plt.show()
    
    def propagate(self):      

        for i in range(len(self.test_particles_list)):
            self.test_particles_list[i] = np.dot(self.transfer_matrix, self.test_particles_list[i])
        
        for i in range(len(self.mismatched_particles_list)):
            self.mismatched_particles_list[i] = np.dot(self.transfer_matrix, self.mismatched_particles_list[i])
        
        for i in range(len(self.matched_particles_list)):
            self.matched_particles_list[i] = np.dot(self.transfer_matrix, self.matched_particles_list[i])

def transfer_matrix(twiss_list):
    mu = twiss_list[3][-1]
    alpha = twiss_list[2][-1]
    beta = twiss_list[1][-1]
    gamma = (1. + alpha**2)/beta 
    
    m11 = np.cos(mu) + np.sin(mu)*alpha
    m12 = beta*np.sin(mu)
    m21 = -gamma*np.sin(mu)
    m22 = np.cos(mu) - np.sin(mu)*alpha
    
    return np.array([[m11,m12],[m21,m22]])
