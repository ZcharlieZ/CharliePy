#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:11:57 2018

@author: charlie
"""
import numpy as np
from scipy.interpolate import interp1d

class Catalogue:
    
    def __init__ (self, xx_part, yy_part, zz_part, mm_part, id_part):
        
        self.xx_part   = xx_part
        
        self.yy_part   = yy_part
        
        self.zz_part   = zz_part
        
        self.mm_part   = mm_part
        
        self.id_part   = id_part
    
    #def __repr__(self):
    #   return repr((self.xx_part, self.yy_part, self.zz_part, self.mm_part, self.id_part))    
        
    def distance(self, xc=0, yc=0, zc=0, a_s=1, b_s=1, c_s=1):
        
        return np.sqrt(((self.xx_part-xc)/a_s)**2+((self.yy_part-yc)/b_s)**2+((self.zz_part-zc)/c_s)**2)
    
    
    
    def centre_of_mass(self, refine=True, stop_number_particle = 1000):
        
        mm_tot = np.sum(self.mm_part)
        
        if refine is not True:
            
            com = (np.sum(self.xx_part*self.mm_part)/mm_tot, np.sum(self.yy_part*self.mm_part)/mm_tot, np.sum(self.zz_part*self.mm_part)/mm_tot)
        
        if refine is True:
            
            xx_temp = self.xx_part
            
            yy_temp = self.yy_part
            
            zz_temp = self.zz_part
            
            mm_temp = self.mm_part
            
            id_temp = self.id_part        
            
            com = (np.sum(self.xx_part*self.mm_part)/mm_tot, np.sum(self.yy_part*self.mm_part)/mm_tot, np.sum(self.zz_part*self.mm_part)/mm_tot)
    
            max_dist_temp = max(self.distance(com[0], com[1], com[2]))
            
            while (len(xx_temp) >= stop_number_particle):
                temp_cat = Catalogue(xx_temp, yy_temp, zz_temp, mm_temp, id_temp)
                
                mm_tot_temp = np.sum(mm_temp)
                
                com  = (np.sum(xx_temp*mm_temp)/mm_tot_temp, np.sum(yy_temp*mm_temp)/mm_tot_temp, np.sum(zz_temp*mm_temp)/mm_tot_temp)

                dist = temp_cat.distance(com[0], com[1], com[2])
    
                xx_temp = np.array(xx_temp[dist<=max_dist_temp])
                
                yy_temp = np.array(yy_temp[dist<=max_dist_temp])
                
                zz_temp = np.array(zz_temp[dist<=max_dist_temp])
                
                mm_temp = np.array(mm_temp[dist<=max_dist_temp])
                
                id_temp = np.array(id_temp[dist<=max_dist_temp])
                
                max_dist_temp = max_dist_temp-0.1*max_dist_temp

        return com
    
    def moment_of_inertia_tensor(self):

        com = self.centre_of_mass(refine=True)
                
        xx_com_shifted = self.xx_part #+ com[0]
        yy_com_shifted = self.yy_part #+ com[1]
        zz_com_shifted = self.zz_part #+ com[2]

        i_xx =  np.sum(xx_com_shifted * xx_com_shifted *self.mm_part)
        i_yy =  np.sum(yy_com_shifted * yy_com_shifted *self.mm_part)
        i_zz =  np.sum(zz_com_shifted * zz_com_shifted *self.mm_part)        
        i_xy =  np.sum(xx_com_shifted * yy_com_shifted * self.mm_part) # I_xy = I_yx
        i_yz =  np.sum(yy_com_shifted * zz_com_shifted * self.mm_part) # I_yz = I_zy
        i_xz =  np.sum(xx_com_shifted * zz_com_shifted * self.mm_part) # I_xz = I_zx

        inertia=np.array([i_xx,i_yy,i_zz])
        eigenvalues, eigenvectors = np.linalg.eigh([[i_xx, i_xy, i_xz], [i_xy, i_yy, i_yz], [i_xz, i_yz, i_zz]])

#        idx = np.argsort(inertia)[::-1]   
#        eigenvalues = eigenvalues[idx]
#        eigenvectors = eigenvectors[:,idx]

        sort_indices = np.argsort(eigenvalues)[::-1]
        c1_v1, c2_v1, c3_v1 = eigenvectors[:, sort_indices[0]]
        c1_v2, c2_v2, c3_v2 = eigenvectors[:, sort_indices[1]]
        
        semiaxis = np.sqrt(eigenvalues)
        print(semiaxis)
        return semiaxis, eigenvalues, eigenvectors
    
    def rotation(self, coords1, coords2):
        
        # WARNING - WORK IN PROGRESS FUNCTION


        coords = np.vstack([coords1, coords2])
        
        cov = np.cov(coords)
        evals, evecs = np.linalg.eig(cov)
        print(evals)
        sort_indices = np.argsort(evals)[::-1]
        c1_v1, c2_v1 = evecs[:, sort_indices[0]]  # Eigenvector with largest eigenvalue
        c1_v2, c2_v2 = evecs[:, sort_indices[1]]

        theta = np.tanh((c1_v1)/(c2_v1))  
        rotation_mat = np.matrix([[np.cos(theta), -np.sin(theta)],\
                                  [np.sin(theta),  np.cos(theta)]])
        transformed_mat = rotation_mat * coords
        
        c1_transformed, c2_transformed = transformed_mat.A
                
        return c1_transformed, c2_transformed, theta, c1_v1, c2_v1, c1_v2, c2_v2


    def sort_for_id(self):        
        
        data = []
        
        for i in range(len(self.id_part)):
            
            data.append(Catalogue(self.xx_part[i], self.yy_part[i], self.zz_part[i], self.mm_part[i], self.id_part[i])) 
            
        data=sorted(data, key=lambda dd: dd.id_part)
        
        return data 
    
    

    def catalogue_properties(self):
        
        print("Particle number : ", len(self.xx_part)    )
        
        print("Minimum distance: ", min(self.distance()) )
        
        print("Maximum distance: ", max(self.distance()) )
        
        print("Centre of mass  : ", self.centre_of_mass(refine=True))
        
        

    def add_a_coloumn(self, number=0):
        
        new_column = np.zeros(len(self.xx_part))+number
        
        return new_column
    
    
    
    def alpha_imf_gradient (dist_1, dist_2, ratio_imf=1, n_bin=200, n_smooth=10):
    
        min_r = min(min(dist_1), min(dist_2))
    
        max_r = max(max(dist_1), max(dist_2))
        
        r_bin=np.logspace(np.log10(min_r), np.log10(max_r), n_bin)
        
        counts1, bins1 = np.histogram(dist_1, r_bin)
        
        counts2, bins2 = np.histogram(dist_2, r_bin)

        n=n_smooth

        temp_counts1 = [np.sum(counts1[(i-n):(i+n+1)])/len(counts1[(i-n):(i+n+1)]) if (i-n) >= 0 else np.sum(counts1[:(i+n+1)])/len(counts1[:(i+n+1)]) for i in range(0,len(counts1))]
        
        temp_counts2 = [np.sum(counts2[(i-n):(i+n+1)])/len(counts2[(i-n):(i+n+1)]) if (i-n) >= 0 else np.sum(counts2[:(i+n+1)])/len(counts2[:(i+n+1)]) for i in range(0,len(counts2))]
    
        counts1 = np.array(temp_counts1)
    
        counts2 = np.array(temp_counts2)    
        
        counts = np.stack((counts1, counts2), axis=-1)
    
        bins = (bins1[1:]+bins1[:-1])/2
        
        counts = np.array([(counts[i][0] + ratio_imf*counts[i][1])/(counts[i][0] + counts[i][1]) for i in range(len(counts))])

        return counts, bins, counts1, counts2
    
    

    def mass_profile(self, n_bin=200, refine=True):
        
        com = self.centre_of_mass(refine=refine)
        
        dist = self.distance(com[0], com[1], com[2])
    
        r_bin = np.logspace(np.log10(min(dist)), np.log10(max(dist)), n_bin)
        
        total_mass = np.sum(self.mm_part)
        
        r_bin = (r_bin[1:]+r_bin[:-1])/2

        # dist_part_w_r = np.array([dist[dist<r_bin[i]] for i in range(len(r_bin))])
        
        mass_part_w_r = np.array([self.mm_part[dist<r_bin[i]] for i in range(len(r_bin))])
        
        mass_w_r = np.array([np.sum(mass_part_w_r[i]) for i in range(len(mass_part_w_r))])
        
        mass_normed = mass_w_r/total_mass
        
        return mass_normed, r_bin
    
    
    
    def density_profile(self, n_bin=200, n_smooth=0, refine=True):
        
        com = self.centre_of_mass(refine=refine)
        
        dist = Catalogue.distance(self, com[0], com[1], com[2])
        
        r_bin = np.logspace(np.log10(min(dist)), np.log10(max(dist)), n_bin)
        
        r_bin = (r_bin[1:]+r_bin[:-1])/2

        counter_0=np.sum([dist<r_bin[0]][0])
        
        density_0=counter_0*self.mm_part[0]/(4./3.*np.pi*(r_bin[0]**3))

        # counter_b_r_r = np.array([np.sum((dist>=r_bin[i-1]) & (dist<r_bin[i])) for i in range(1, len(r_bin))])
        
        mass_part_b_r_r = np.array([self.mm_part[(dist>=r_bin[i-1]) & (dist<r_bin[i])] for i in range(1, len(r_bin))])
        
        mass_b_r_r = np.array([np.sum(mass_part_b_r_r[i]) for i in range(len(mass_part_b_r_r))])
        
        density=np.array([mass_b_r_r[i]/(4./3.*np.pi*(r_bin[i]**3 - r_bin[i-1]**3))  for i in range(len(mass_b_r_r))])
        
        density=np.insert(density, 0, density_0)
                
        n=n_smooth

        temp_density = [np.sum(density[(i-n):(i+n+1)])/len(density[(i-n):(i+n+1)]) if (i-n) >= 0 else np.sum(density[:(i+n+1)])/len(density[:(i+n+1)]) for i in range(0,len(density))]

        density = temp_density

        return density, r_bin
    
    

    def radius_given_mass(self, percentage=50, refine=True):
        
        m_profile, r_profile = self.mass_profile(refine=refine)
        
        function = interp1d (m_profile, r_profile)
        
        return function(percentage/100)