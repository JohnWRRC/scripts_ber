# Loading libraries
import sys, os
import numpy as np
from scipy.stats import lognorm, gamma
import time

def create_environment(env_dim = 10000, plant_sp = 3, individuals = 1500, dist_corr = 200.0, rho = 0.0, stage = 0,
                       file_name = ""):
                       
                       
    '''
    This function populates a continuous landscape of dimension env_dim (m) with "individuals" 
    plants of "plant_sp" different species. dist_corr and rho sets the autocorralation in the position
    of plants.
    
    We still need:
    - aggregate species in clusters?
    '''
    
    
    #vars
    cdef float pl_sp,x,y,names,order,prop,spid,index_ind,propacum
    cdef int s,nind
    cdef c_array.array coords = np.array([[], []], dtype = float)
    
    
    # Plant agents
    class Plants:
        #vars
        def __init__( self):
            cdef float idd,sp
            cdef c_array.array coords = np.array([[], []], dtype = float)
            cdef c_array.array fruits = np.array([], dtype = int)

        idd = np.arange(1, individuals+1)
        
        #coords = np.array([[], []], dtype = float)
        sp = np.zeros((individuals,), dtype=int)
        
    # instanciando a classe plants    
    plants = Plants()
    
    
    pl_sp = plant_sp
    
    #np.random.seed(3)
    
    # Plant location
    x = env_dim*np.random.random(individuals)
    y = env_dim*np.random.random(individuals)
    plants.coords = np.append(x, y).reshape(2, individuals).transpose()
    
    # Number of fruits
    
    plants.fruits = gamma.rvs(2, scale = 8, size = individuals).astype(int)
    np.putmask(plants.fruits, plants.fruits > 100, 100)
    
    
    # Defining plant species
    if Parms.p_abund_field:
        names, prop, mat = read_abund(file_name, stage)
        order = np.argsort(prop)[::-1]
        names = names[order]
        prop = np.sort(prop)[::-1]
    else:
        names = np.arange(1, pl_sp+1).astype('S10')
        spid = np.arange(1, pl_sp+1) # species id
        s = 2 # parameter of the lognormal distribution
        prop = lognorm.pdf(spid, s, loc=0, scale=1) # proportions come from a discretized log normal distribution, with parameters (mulog = 0, sdlog = 2)

    prop = prop/sum(prop) # proportion of individuals of each species
    nind = np.rint(prop*individuals).astype(int) # number of individuals of each species
    # removing sps not present
    if Parms.p_remove_absent:    
        names = names[np.where(nind != 0)]
        nind = nind[np.where(nind != 0)]
        pl_sp = len(nind)
        
    while nind.sum() < individuals:
        x = np.random.choice(np.arange(pl_sp))
        nind[x] = nind[x] + 1
    while nind.sum() > individuals:
        x = np.random.choice(np.arange(pl_sp))
        nind[x] = nind[x] - 1
    prop = nind.astype(float)/nind.sum() # proportion of individuals of each species
    propacum = prop.cumsum(0) # Cumulative probability for each species
    
    
    #vars
    cdef float nrand,index_sp,prop_aux
   
    
    # First plant
    x = np.where( plants.sp == 0 )[0]
    index_ind = np.random.choice(x)
    nrand = np.random.uniform()
    index_sp = np.amin(np.where( nrand < propacum ))
    
    plants.sp[index_ind] = index_sp+1
    
    # Refreshing
    nind[index_sp] = nind[index_sp]-1
    prop = nind.astype(float)/nind.sum()
    
    if prop[index_sp] > 0.0:
        prop_aux = np.delete(prop, index_sp)
        prop_aux = np.array(map(lambda x: x * (1 - rho), prop_aux))
        prop_aux = np.insert(prop_aux, index_sp, prop[index_sp] + (1 - prop[index_sp]) * rho)
        propacum = prop_aux.cumsum(0)
    else:
        propacum = prop.cumsum(0)
    
    # Other plants
    #while np.any(plants.sp == 0):
    cdef c_array.array dists = np.array([[], []], dtype = float)
    cdef float dist_index
    while nind.sum() > 0:
        
        dists = np.array(map(distance, np.repeat(plants.coords[index_ind].reshape((1,2)), individuals, axis=0), plants.coords))
        
        dist_index = np.where( (dists < dist_corr) * (plants.sp == 0) )[0]
        if dist_index.size != 0:
            index_ind = np.random.choice(dist_index)
        else:
            #dists_sort = np.sort(dists) 
            #for i in np.arange(len(dists)):
            #    indice = np.where( dists_sort[i] == dists )
            #    if plants.sp[indice] == 0:
            #        index_ind = indice
            #        break;
            x = np.where( plants.sp == 0 )[0]
            index_ind = np.random.choice(x)
            
        nrand = np.random.uniform()
        index_sp = np.amin(np.where( nrand < propacum ))
        
        plants.sp[index_ind] = index_sp+1
        
        # Refreshing
        nind[index_sp] = nind[index_sp]-1
        prop = nind.astype(float)/nind.sum()
        print nind.sum()
        
        if prop[index_sp] > 0.0:
            prop_aux = np.delete(prop, index_sp)
            prop_aux = np.array(map(lambda x: x * (1 - rho), prop_aux))
            prop_aux = np.insert(prop_aux, index_sp, prop[index_sp] + (1 - prop[index_sp]) * rho)
            propacum = prop_aux.cumsum(0)
        else:
            propacum = prop.cumsum(0)
            
    return plants.idd, plants.sp, plants.fruits, plants.coords, names