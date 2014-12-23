
# a simple set of functions to compute monte carlo simulations of a lattice gas using the metropolis algorithm
# William Gilpin 2014


import numpy
import matplotlib
import scipy

def neighbor_engy(lat, ii, jj):
    (xspan, yspan) = lat.shape
    engy = 0.0
    if ii == 0:
        engy += -.5*lat[ii, jj]*lat[xspan-1, jj]
    else:
        engy += -.5*lat[ii, jj]*lat[ii-1, jj]


    if jj == 0:
        engy += -.5*lat[ii, jj]*lat[ii, yspan-1]
    else:
        engy += -.5*lat[ii, jj]*(lat[ii, jj-1])

    if ii == (xspan-1):
        engy += -.5*lat[ii, jj]*(lat[0, jj])
    else:
        engy += -.5*lat[ii, jj]*(lat[ii+1, jj])

    if jj == (yspan-1):
        engy += -.5*lat[ii, jj]*(lat[ii, 0])
    else:
        engy += -.5*lat[ii, jj]*(lat[ii, jj+1])
    return engy

def lattice_energy(lat):
    (xspan, yspan) = lat.shape
    engy = 0.0
    for ii in xrange(xspan):
        for jj in xrange(yspan):
            engy += neighbor_engy(lat, ii, jj)
    return engy

# only move to nearest neighbor
def update_lattice(lat0, T):
    lat = copy(lat0)
    (xspan, yspan) = lat.shape

    engy0 = lattice_energy(lat)

    (i_pt, j_pt) = (floor(xspan*rand(1))[0], floor(xspan*rand(1))[0])
    swt = rand(1)
    if (swt < .5):
        i_pt2 = mod(i_pt+1, xspan-1)
        j_pt2 = j_pt   
    else:
        i_pt2 = i_pt
        j_pt2 = mod(j_pt+1, yspan-1)
    
    if lat[i_pt, j_pt] == (1 - lat[i_pt2, j_pt2]):
        lat[i_pt, j_pt] = 1 - lat[i_pt, j_pt]
        lat[i_pt2, j_pt2] = 1 - lat[i_pt2, j_pt2]
    elif lat[i_pt, j_pt] == lat[i_pt2, j_pt2]:
        pass
    
    engy1 = lattice_energy(lat)

    del_engy = engy1 - engy0
    
    if del_engy <= 0:
        pass
    else:
        die = rand(1)
        if die <= exp(-del_engy/T):
            pass
        else:
            lat[i_pt, j_pt] = 1 - lat[i_pt, j_pt]
            lat[i_pt2, j_pt2] = 1 - lat[i_pt2, j_pt2]
    return lat

def mc_anneal(lat0, T, nsteps=10000):
    lat = copy(lat0)
    for ii in xrange(nsteps):
        lat = copy(update_lattice(lat, T))
    lat2 = lat
    return lat2

# number of points in the lattice
num_pts = 40;

# different temperatures at which to perform simulation
temps = linspace(.01, 5, num_pts)

engys = zeros(num_pts)
engy_sq = zeros(num_pts)
heat_cap = zeros(num_pts)

# Run metropolis
for ii in xrange(len(temps)):
    temp = temps[ii]
    
    # number of runs at that temp to average over
    N = 20
    engy_temp = list()
    engysq_temp = list()
    part_func_temp = list()
    for jj in xrange(N): 
        lattice = floor(2*rand(10,10))
        lattice2 = mc_anneal(lattice, temp)
        lat_engy = lattice_energy(lattice2)
        engy_temp.append( lat_engy*exp(-lat_engy/temp))
        engysq_temp.append( (lat_engy**2)*exp(-lat_engy/temp)  )
        part_func_temp.append( exp(-lat_engy/temp) )
                           
    part_func = sum(part_func_temp)
    engys[ii] = sum(engy_temp)/part_func
    engy_sq[ii] = sum(engysq_temp)/part_func

    heat_cap[ii] = (engy_sq[ii] - engys[ii]**2)/(temp**2)
    
    
    print (str(ii) + ", temperature: " + str(temp))
      
# A discontinuity appears at the critical temperature
plot(temps, heat_cap,'.r',markersize=10)
xlabel("Temperature")
ylabel("Heat Capacity")
figure()
semilogy(temps, heat_cap,'.r',markersize=10)
xlabel("Temperature")
ylabel("Heat Capacity")



