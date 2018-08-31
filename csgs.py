import numpy
import time





#For all five functions, x,y,z are one-dimensional arrays of floats with coordinates in micrometers, f is a float with the
#equivalent focal of the system in mm, d is a float with the pixel size in micrometers, lam is the wavelength used in
#micrometers, res is the lateral resolution of the SLM.The software assumes the use of a circular area within a square slm:
#additional zero padding may be needed to operate rectangular SLMs

# Random superposition algorithm: Fastest available algorithm, produces low quality holograms

def rs(x,y,z,f,d,lam,res):
    t=time.clock()

    #creation of a list of the SLM pixels contained in the pupil
    slm_xcoord,slm_ycoord=numpy.meshgrid(numpy.linspace(-1.0,1.0,res),numpy.linspace(-1.0,1.0,res))
    pup_coords=numpy.where(slm_xcoord**2+slm_ycoord**2<1.0)

    #array containing the phase of the field at each created spot
    pists=numpy.random.random(x.shape[0])*2*numpy.pi

    #conversion of the coordinates arrays in microns
    slm_xcoord=slm_xcoord*d*float(res)/2.0
    slm_ycoord=slm_ycoord*d*float(res)/2.0
    
    #computation of the phase patterns generating each single spot independently
    slm_p_phase=numpy.zeros((x.shape[0],pup_coords[0].shape[0]))
    for i in range(x.shape[0]):
        slm_p_phase[i,:]=2.0*numpy.pi/(lam*(f*10.0**3))*(x[i]*slm_xcoord[pup_coords]+y[i]*slm_ycoord[pup_coords])+(numpy.pi*z[i])/(lam*(f*10.0**3)**2)*(slm_xcoord[pup_coords]**2+slm_ycoord[pup_coords]**2)


    #creation of the hologram, as superposition of all the phase patterns with random pistons
    slm_total_field=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_p_phase+pists[:,None])),axis=0)
    slm_total_phase=numpy.angle(slm_total_field)

    t=time.clock()-t

    #evaluation of the algorithm performance, calculating the expected intensities of all spots

    spot_fields=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_total_phase[None,:]-slm_p_phase)),axis=1)
    ints=numpy.abs(spot_fields)**2


    #reshaping of the hologram in a square array

    out=numpy.zeros((res,res))
    out[pup_coords]=slm_total_phase

    #the function returns the hologram, and a list with efficiency, uniformity and variance of the spots, and hologram computation time
    
    return out,[numpy.sum(ints),1-(numpy.max(ints)-numpy.min(ints))/(numpy.max(ints)+numpy.min(ints)),numpy.sqrt(numpy.var(ints))/numpy.mean(ints),t]




# Standard GS algorithm: Slow, high efficiency holograms, better uniformity than RS. The parameter "iters" is the number of GS iterations to
# perform

def gs(x,y,z,f,d,lam,res,iters):
    t=time.clock()

    #creation of a list of the SLM pixels contained in the pupil
    slm_xcoord,slm_ycoord=numpy.meshgrid(numpy.linspace(-1.0,1.0,res),numpy.linspace(-1.0,1.0,res))
    pup_coords=numpy.where(slm_xcoord**2+slm_ycoord**2<1.0)
    
    
    #array containing the phase of the field at each created spot
    pists=numpy.random.random(x.shape[0])*2*numpy.pi

    #conversion of the coordinates arrays in microns
    slm_xcoord=slm_xcoord*d*float(res)/2.0
    slm_ycoord=slm_ycoord*d*float(res)/2.0
    
    #computation of the phase patterns generating each single spot independently
    slm_p_phase=numpy.zeros((x.shape[0],pup_coords[0].shape[0]))
    for i in range(x.shape[0]):
        slm_p_phase[i,:]=2.0*numpy.pi/(lam*(f*10.0**3))*(x[i]*slm_xcoord[pup_coords]+y[i]*slm_ycoord[pup_coords])+(numpy.pi*z[i])/(lam*(f*10.0**3)**2)*(slm_xcoord[pup_coords]**2+slm_ycoord[pup_coords]**2)


    #main GS loop
    for n in range(iters):
        #creation of the hologram, as superposition of all the phase patterns with random pistons
        slm_total_field=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_p_phase+pists[:,None])),axis=0)
        slm_total_phase=numpy.angle(slm_total_field)


        #Update of the phases at all spots locations. The intensities are evaluated too for performance
        #estimation of the algorithm
        spot_fields=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_total_phase[None,:]-slm_p_phase)),axis=1)
        pists=numpy.angle(spot_fields)
        ints=numpy.abs(spot_fields)**2

    t=time.clock()-t
            
    #the function returns the hologram, and a list with efficiency, uniformity and variance of the spots, and hologram computation time
    out=numpy.zeros((res,res))
    out[pup_coords]=slm_total_phase
    return out,[numpy.sum(ints),1-(numpy.max(ints)-numpy.min(ints))/(numpy.max(ints)+numpy.min(ints)),numpy.sqrt(numpy.var(ints))/numpy.mean(ints),t]



# Standard WGS algorithm: Slow, high uniformity holograms, better efficiency than RS. The parameter "iters" is the number of GS iterations to
# perform

def wgs(x,y,z,f,d,lam,res,iters):
    t=time.clock()

    #creation of a list of the SLM pixels contained in the pupil
    slm_xcoord,slm_ycoord=numpy.meshgrid(numpy.linspace(-1.0,1.0,res),numpy.linspace(-1.0,1.0,res))
    pup_coords=numpy.where(slm_xcoord**2+slm_ycoord**2<1.0)
    
    #initialization of the weights, all with equal value
    weights=numpy.ones(x.shape[0])/float(x.shape[0])
    
    #array containing the phase of the field at each created spot
    pists=numpy.random.random(x.shape[0])*2*numpy.pi

    #conversion of the coordinates arrays in microns
    slm_xcoord=slm_xcoord*d*float(res)/2.0
    slm_ycoord=slm_ycoord*d*float(res)/2.0
    
    #computation of the phase patterns generating each single spot independently
    slm_p_phase=numpy.zeros((x.shape[0],pup_coords[0].shape[0]))
    for i in range(x.shape[0]):
        slm_p_phase[i,:]=2.0*numpy.pi/(lam*(f*10.0**3))*(x[i]*slm_xcoord[pup_coords]+y[i]*slm_ycoord[pup_coords])+(numpy.pi*z[i])/(lam*(f*10.0**3)**2)*(slm_xcoord[pup_coords]**2+slm_ycoord[pup_coords]**2)


    #main GS loop
    for n in range(iters):
        #Creation of the hologram, as superposition of all the phase patterns with random pistons and weighted intensity
        slm_total_field=numpy.sum(weights[:,None]/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_p_phase+pists[:,None])),axis=0)
        slm_total_phase=numpy.angle(slm_total_field)


        #Update of the phases at all spots locations. The intensities are evaluated too for performance
        #estimation of the algorithm
        spot_fields=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_total_phase[None,:]-slm_p_phase)),axis=1)
        pists=numpy.angle(spot_fields)
        ints=numpy.abs(spot_fields)**2
        #update and renormalization of the weights. Renormalization is required to avoid computation precision errors in the final performance for
        #large numbers of iterations
        weights=weights*(numpy.mean(numpy.sqrt(ints))/numpy.sqrt(ints))
        weights=weights/numpy.sum(weights)


    t=time.clock()-t
            
    #the function returns the hologram, and a list with efficiency, uniformity and variance of the spots, and hologram computation time
    out=numpy.zeros((res,res))
    out[pup_coords]=slm_total_phase
    return out,[numpy.sum(ints),1-(numpy.max(ints)-numpy.min(ints))/(numpy.max(ints)+numpy.min(ints)),numpy.sqrt(numpy.var(ints))/numpy.mean(ints),t]



# Compressive sensing GS:Fast, similar performances to GS. The parameter "sub" is a float greater than 0.0 and up to 1.0, indicating the
# Subsampling of the pupil. For sub=1.0, CSGS is equivalent to GS. Smaller values of sub increase the speed, with a maximum speed
# equal to the speed of RS. If sub is too small, performance may be affected (as a rule of thumb res^2*sub should be bigger than the
# number of spots)

def csgs(x,y,z,f,d,lam,res,iters,sub):
    t=time.clock()
    #creation of a list of the SLM pixels contained in the pupil
    slm_xcoord,slm_ycoord=numpy.meshgrid(numpy.linspace(-1.0,1.0,res),numpy.linspace(-1.0,1.0,res))
    pup_coords=numpy.where(slm_xcoord**2+slm_ycoord**2<1.0)

    #creation of a list of the indexes of pup_coords, shuffled in random order
    coordslist=numpy.asarray(range(pup_coords[0].shape[0]))
    numpy.random.shuffle(coordslist)
    
    #array containing the phase of the field at each created spot
    pists=numpy.random.random(x.shape[0])*2*numpy.pi

    
    #conversion of the coordinates arrays in microns
    slm_xcoord=slm_xcoord*d*float(res)/2.0
    slm_ycoord=slm_ycoord*d*float(res)/2.0
    
    #computation of the phase patterns generating each single spot independently
    slm_p_phase=numpy.zeros((x.shape[0],pup_coords[0].shape[0]))
    for i in range(x.shape[0]):
        slm_p_phase[i,:]=2.0*numpy.pi/(lam*(f*10.0**3))*(x[i]*slm_xcoord[pup_coords]+y[i]*slm_ycoord[pup_coords])+(numpy.pi*z[i])/(lam*(f*10.0**3)**2)*(slm_xcoord[pup_coords]**2+slm_ycoord[pup_coords]**2)

    #main GS loop
    for n in range(iters):
        #a new set of random points is chosen on the SLM
        coordslist=numpy.roll(coordslist,int(coordslist.shape[0]*sub))
        coordslist_sparse=coordslist[:int(coordslist.shape[0]*sub)]

        #Creation of the hologram, as superposition of all the phase patterns with the estimated phases and weighted intensity
        slm_total_field=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_p_phase[:,coordslist_sparse]+pists[:,None])),axis=0)
        slm_total_phase=numpy.angle(slm_total_field)

        #Update of the phases at all spots locations. The intensities are evaluated too for performance
        #estimation of the algorithm
        spot_fields=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_total_phase[None,:]-slm_p_phase[:,coordslist_sparse])),axis=1)
        pists=numpy.angle(spot_fields)
        ints=numpy.abs(spot_fields)**2


    #The output holograms is estimated with the full SLM resolution
    slm_total_field=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_p_phase+pists[:,None])),axis=0)
    slm_total_phase=numpy.angle(slm_total_field)

    t=time.clock()-t


    #evaluation of the algorithm performance, calculating the expected intensities of all spots
    spot_fields=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_total_phase[None,:]-slm_p_phase)),axis=1)
    ints=numpy.abs(spot_fields)**2
            
    #the function returns the hologram, and a list with efficiency, uniformity and variance of the spots, and hologram computation time
    out=numpy.zeros((res,res))
    out[pup_coords]=slm_total_phase
    return out,[numpy.sum(ints),1-(numpy.max(ints)-numpy.min(ints))/(numpy.max(ints)+numpy.min(ints)),numpy.sqrt(numpy.var(ints))/numpy.mean(ints),t]

# Weighted compressive sensing GS:Fast, efficiency and uniformity between GS and WGS.
# The parameter "sub" is a float greater than 0.0 and up to 1.0, indicating the Subsampling of the pupil. For sub=1.0, CSGS is equivalent to GS.
# Smaller values of sub increase the speed, with a maximum speed equal to half the speed of RS. If sub is too small, performance may be affected
# (as a rule of thumb (res^2)*sub should be at least twice the number of spots)

def wcsgs(x,y,z,f,d,lam,res,iters,sub):
    t=time.clock()
    #creation of a list of the SLM pixels contained in the pupil
    slm_xcoord,slm_ycoord=numpy.meshgrid(numpy.linspace(-1.0,1.0,res),numpy.linspace(-1.0,1.0,res))
    pup_coords=numpy.where(slm_xcoord**2+slm_ycoord**2<1.0)

    #creation of a list of the indexes of pup_coords, shuffled in random order
    coordslist=numpy.asarray(range(pup_coords[0].shape[0]))
    numpy.random.shuffle(coordslist)
    
    #array containing the phase of the field at each created spot
    pists=numpy.random.random(x.shape[0])*2*numpy.pi

    
    #conversion of the coordinates arrays in microns
    slm_xcoord=slm_xcoord*d*float(res)/2.0
    slm_ycoord=slm_ycoord*d*float(res)/2.0
    
    #computation of the phase patterns generating each single spot independently
    slm_p_phase=numpy.zeros((x.shape[0],pup_coords[0].shape[0]))
    for i in range(x.shape[0]):
        slm_p_phase[i,:]=2.0*numpy.pi/(lam*(f*10.0**3))*(x[i]*slm_xcoord[pup_coords]+y[i]*slm_ycoord[pup_coords])+(numpy.pi*z[i])/(lam*(f*10.0**3)**2)*(slm_xcoord[pup_coords]**2+slm_ycoord[pup_coords]**2)

    #main GS loop
    for n in range(iters-1):
        #a new set of random points is chosen on the SLM
        coordslist=numpy.roll(coordslist,int(coordslist.shape[0]*sub))
        coordslist_sparse=coordslist[:int(coordslist.shape[0]*sub)]

        #Creation of the hologram, as superposition of all the phase patterns with the estimated phases
        slm_total_field=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_p_phase[:,coordslist_sparse]+pists[:,None])),axis=0)
        slm_total_phase=numpy.angle(slm_total_field)

        #Update of the phases at all spots locations. The intensities are evaluated too for performance
        #estimation of the algorithm
        spot_fields=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_total_phase[None,:]-slm_p_phase[:,coordslist_sparse])),axis=1)
        pists=numpy.angle(spot_fields)
        ints=numpy.abs(spot_fields)**2


    # an additional single loop of WGS without compression
    slm_total_field=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_p_phase+pists[:,None])),axis=0)
    slm_total_phase=numpy.angle(slm_total_field)

    spot_fields=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_total_phase[None,:]-slm_p_phase)),axis=1)
    pists=numpy.angle(spot_fields)
    ints=numpy.abs(spot_fields)**2

    weights=numpy.ones(x.shape[0])/float(x.shape[0])*(numpy.mean(numpy.sqrt(ints))/numpy.sqrt(ints))
    weights=weights/numpy.sum(weights)

    slm_total_field=numpy.sum(weights[:,None]/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_p_phase+pists[:,None])),axis=0)
    slm_total_phase=numpy.angle(slm_total_field)

    t=time.clock()-t

    #evaluation of the algorithm performance, calculating the expected intensities of all spots
    spot_fields=numpy.sum(1.0/(float(pup_coords[0].shape[0]))*numpy.exp(1j*(slm_total_phase[None,:]-slm_p_phase)),axis=1)
    ints=numpy.abs(spot_fields)**2    
            
    #the function returns the hologram, and a list with efficiency, uniformity and variance of the spots, and hologram computation time
    out=numpy.zeros((res,res))
    out[pup_coords]=slm_total_phase
    return out,[numpy.sum(ints),1-(numpy.max(ints)-numpy.min(ints))/(numpy.max(ints)+numpy.min(ints)),numpy.sqrt(numpy.var(ints))/numpy.mean(ints),t]

if __name__=="__main__":

    x=(numpy.random.random(100)-0.5)*100.0
    y=(numpy.random.random(100)-0.5)*100.0
    z=(numpy.random.random(100)-0.5)*10.0

    performance_pars=["Efficiency : ","Uniformity : ","Variance : ","Computation time (s) : "]

    print "Computing random superposition hologram:"

    phase, performance=rs(x,y,z,20.0,15.0,0.488,512)

    for i in range(4):
        print performance_pars[i],performance[i]

    print "Computing Gerchberg-Saxton hologram:"

    phase, performance=gs(x,y,z,20.0,15.0,0.488,512,30)
    
   
    for i in range(4):
        print performance_pars[i],performance[i]

    print "Computing Weighted Gerchberg-Saxton hologram:"

    phase, performance=wgs(x,y,z,20.0,15.0,0.488,512,30)

    
    for i in range(4):
        print performance_pars[i],performance[i]

    print "Computing Compressive Sensing Gerchberg-Saxton hologram:"

    phase, performance=csgs(x,y,z,20.0,15.0,0.488,512,30,0.05)


    for i in range(4):
        print performance_pars[i],performance[i]

    print "Computing Weighted Compressive Sensing Gerchberg-Saxton hologram:"

    phase, performance=wcsgs(x,y,z,20.0,15.0,0.488,512,30,0.05)

    for i in range(4):
        print performance_pars[i],performance[i]

