
import numpy as np
import os, sys
import finufft
from vtk.util import numpy_support
from vtk import vtkXMLPPolyDataReader

pi = np.pi

"""
Locate source and target data directories

    Inputs:
        nu (float) : volume fraction
        L (int) : box size
        directories (list) : list of parent directories to search
        target (str) : target directory for processed data
    Outputs:
        source (str) : source directory
        target (str) : target directory

"""
def locateData(nu, L, directories, target):

    # Source directory
    source = f"v1.0_nu{nu}_L{L}_A5"

    # Check which parent directory
    for directory in directories:
        if os.path.exists(f"{directory}/{source}"):
            print(f"Found source directory in {directory}")
            source = f"{directory}/{source}"
            break

    # If not found, exit
    if not os.path.exists(f"{source}"):
        print(f"Source directory {source} does not exist!")
        sys.exit(1)

    # Make target directory
    target = f"{target}/nu{nu}_L{L}"
    try: os.makedirs(target)
    except: pass

    return source, target

"""
Returns sorted array of particle positions and orientations

    Inputs:
        source (str) : filename

    Outputs:
        X (N x 3 array) : sorted particle positions
        P (N x 3 array) : sorted particle orientations
"""
def configuration_sorted(source):

    # Load in data
    reader = vtkXMLPPolyDataReader()
    reader.SetFileName(source)
    reader.Update()

    # Get fields
    cell_data = reader.GetOutput().GetCellData()
    field_names = [cell_data.GetArrayName(i) for i in range(cell_data.GetNumberOfArrays())]

    # Get global particle id
    fid = field_names.index("gid")
    globalId = numpy_support.vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(fid))

    # Sort
    sortedId = np.argsort(globalId)
    
    # Get points
    X = numpy_support.vtk_to_numpy(reader.GetOutput().GetPoints().GetData()).astype(np.float64)
    X = 0.5 * (X[1::2] + X[0::2])
    X = X[sortedId]

    # Get orientations
    fid = field_names.index("znorm")
    P = numpy_support.vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(fid))
    P = P[sortedId]

    return X, P

"""
Loads particle configuration from vtk files

    Inputs:
            source (str) : filename
            M (int) : number of quadrature points

    Outputs:
            X (N x 3 array) : center of mass coordinates
            P (N x 3 array) : orientation vector
            Y (N * M x 3 array) : quadrature points
            F (N * M x 3 array) : force distribution
"""
def configuration(source, M):

    # Load in data
    reader = vtkXMLPPolyDataReader()
    reader.SetFileName(source)
    reader.Update()

    # Get points
    X = numpy_support.vtk_to_numpy(reader.GetOutput().GetPoints().GetData()).astype(np.float64)

    # Copy for NUFFT
    Y = X.copy()

    def GetPointData(name):
        return numpy_support.vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(name))

    # Get quadrature points
    sQuad = GetPointData("sQuad")

    # Get quadrature weights
    weightQuad = GetPointData("weightQuad")

    # Get forces along quadrature points
    forceHydro = GetPointData("forceHydro")

    # Compute orientation vector
    P = X[M-1::M, :] - X[0::M, :]
    P /= np.linalg.norm(P, axis=1)[:, None]

    # Compute center of mass (symmetric discretization)
    X = 0.5 * (X[M-1::M, :] + X[0::M, :])

    # NUFFT sources
    F = (forceHydro * weightQuad[:, None]).astype(np.float64)

    return X, P, Y, F

"""
Get particle velocities

    Inputs:
        source (str) : filename

    Outputs:
        U (N x 1 array) : particle velocities
"""
def particleVelocity(source):

    # Load in data
    reader = vtkXMLPPolyDataReader()
    reader.SetFileName(source)
    reader.Update()
    
    # Get fields
    number_of_fields = reader.GetOutput().GetCellData().GetNumberOfArrays()
    field_names = [reader.GetOutput().GetCellData().GetArrayName(i) for i in range(number_of_fields)]

    # Compute particle velocity
    fid = field_names.index("znorm")
    znorm = numpy_support.vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(fid))

    # Compute particle velocity
    fid = field_names.index("vel")
    U = numpy_support.vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(fid))

    # Compute particle velocity
    fid = field_names.index("velCollision")
    Uc = numpy_support.vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(fid))

    # Subtract off intrinsic swimming speed
    U = U - znorm

    # Compute hydro velocity
    Uf = np.linalg.norm(U - Uc, axis=1)

    # Compute collision velocity
    Uc = np.linalg.norm(Uc, axis=1)

    # Compute full velocity
    U = np.linalg.norm(U, axis=1)
    
    return np.array([U, Uf, Uc])

"""
Computes the radially averaged autocorrelation function for order parameters and velocity

    Inputs:
        X (N x 3 array) : center of mass coordinates
        P (N x 3 array) : orientation vector
        Y (N * M x 3 array) : quadrature points
        F (N * M x 3 array) : force distribution
        L (int) : box size

    Outputs : 
        r : radial coordinate
        c : concentration correlation function
        n : polarity correlation function
        Q : nematic correlation function
        u : velocity correlation function
"""
def correlation(X, P, Y, F, L):

    # Convert to [0, 2 * pi] for nufft
    k0 = 2 * pi / L #fundamental mode

    x0 = k0 * X[:, 0]
    x1 = k0 * X[:, 1]
    x2 = k0 * X[:, 2]

    y0 = k0 * Y[:, 0]
    y1 = k0 * Y[:, 1]
    y2 = k0 * Y[:, 2]

    # Number of modes
    Nmodes = L // 2
    shp = (Nmodes, Nmodes, Nmodes)

    # Number of particles
    Nparticles = X.shape[0] 

    # Initialize nufft for conformation variables
    plan = finufft.Plan(1, shp, n_trans=1, modeord=1, dtype=np.complex128)
    plan.setpts(x0, x1, x2)

    # Define nufft function
    nufft = lambda f: plan.execute(f.astype(np.complex128)) / Nparticles

    # Initialize nufft for force
    planForce = finufft.Plan(1, shp, n_trans=1, modeord=1, dtype=np.complex128)
    planForce.setpts(y0, y1, y2)
    
    # Construct Fourier modes
    kx = np.fft.fftfreq(Nmodes, d=1.0/(Nmodes * k0))
    kx, ky, kz = np.meshgrid(kx, kx, kx, indexing='ij')

    # Construct wave vectors
    k = np.sqrt(kx**2 + ky**2 + kz**2)

    # Avoid division by zero
    k[0,0,0] = 1

    # Unit wave vectors
    khat = np.stack((kx/k, ky/k, kz/k), axis=-1)

    # Preallocate
    u_h = np.zeros((Nmodes, Nmodes, Nmodes, 3), dtype=np.complex128)

    # Compute Fourier transform of velocity
    for m in range(3):
        f_h = planForce.execute(F[:, m].astype(np.complex128)) / L**3
        for n in range(3):
            u_h[...,n] += ((np.float64(n == m) - khat[...,n] * khat[...,m]) / k ** 2) * f_h

    # Set zero mean
    u_h[0,0,0,:] = 0

    # Concentration
    c_h = nufft(np.ones(len(x0)))
    Pc = np.abs(c_h)**2

    # Preallocate
    Pu, Pn, PQ = (np.zeros(shp) for _ in range(3))

    # Velocity, polarity, and nematic tensor
    for n in range(3):

        Pu += np.abs(u_h[..., n])**2 #velocity
        Pn += np.abs(nufft(P[:, n]))**2 #polarity

        for m in range(3):
            PQ += np.abs(nufft(P[:, n] * P[:, m] - (n == m) / 3))**2 #nematic

    # Compute correlation functions from power spectrum
    c,_ = power2corr(Pc, L)
    n,_ = power2corr(Pn, L)
    Q,_ = power2corr(PQ, L)
    u,r = power2corr(Pu, L)
    
    # Output data
    return np.array([r, c, n, Q, u])

"""
Computes autocorrelation function, normalized by domain size, from power spectrum.

    Inputs:
        P_h : power spectrum
        L : box size
    Outputs : 
        C : radially averaged real-space correlation function
        r : radial coordinate
"""
def power2corr(P_h, L):

    # Get number of Fourier modes
    Nmodes = P_h.shape[0]

    # Compute real-space correlation function via inverse FFT
    C = np.fft.ifftn(P_h) * (Nmodes ** 3)
    C = np.real(np.fft.ifftshift(C)) #shift zero to center

    # Create grids of indices and bin radius to nearest integer
    x, y, z = np.indices((Nmodes, Nmodes, Nmodes)) - (Nmodes // 2)
    r = np.round(np.sqrt(x**2 + y**2 + z**2)).astype(np.int64)

    # Flatten arrays
    r, C = r.ravel(), C.ravel()

    # Bin correlation function
    C = np.bincount(r, weights=C) / np.bincount(r)

    # Compute radial distance
    r = np.arange(len(C)) * (L / Nmodes)

    # Return data
    return C, r

"""
Get particle velocities

    Inputs:
        source (str) : filename

    Outputs:
        U (N x 1 array) : particle velocities
"""
def particleVelocity(source):

    # Load in data
    reader = vtkXMLPPolyDataReader()
    reader.SetFileName(source)
    reader.Update()
    
    # Get fields
    number_of_fields = reader.GetOutput().GetCellData().GetNumberOfArrays()
    field_names = [reader.GetOutput().GetCellData().GetArrayName(i) for i in range(number_of_fields)]

    # Compute particle velocity
    fid = field_names.index("znorm")
    znorm = numpy_support.vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(fid))

    # Compute particle velocity
    fid = field_names.index("vel")
    U = numpy_support.vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(fid))

    # Compute particle velocity
    fid = field_names.index("velCollision")
    Uc = numpy_support.vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(fid))

    # Subtract off intrinsic swimming speed
    U = U - znorm

    # Compute hydro velocity
    Uf = np.linalg.norm(U - Uc, axis=1)

    # Compute collision velocity
    Uc = np.linalg.norm(Uc, axis=1)

    # Compute full velocity
    U = np.linalg.norm(U, axis=1)
    
    return np.array([U, Uf, Uc])