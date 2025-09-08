
# Import modules
import numpy as np
import glob, sys, yaml
from natsort import natsorted
import sbt

from config import directories, target

# Get box size and volume fraction from command line input
try: nu, L = float(sys.argv[1]), int(sys.argv[2]) 
except: 
    print("Invalid input. Please provide volume fraction (nu) and box size (L).")
    sys.exit(1)
	
# Locate filepath
source, target = sbt.locateData(nu, L, directories, f"{target}/msd")

# Get list of files and sort
filenames = glob.glob(f"{source}/result/result*/Sylinder_*.pvtp")
filenames = natsorted(filenames) #sorted list of files
Nfiles = len(filenames) #number of files

# Range of files to process
n0 = Nfiles - 100 if Nfiles > 100 else 0
n1 = Nfiles

# Initial configuration
X0, P0 = sbt.configuration_sorted(filenames[n0])

# Initialize
dX, dP = np.zeros(Nfiles), np.zeros((Nfiles, 3, 3))

# Loop over all files
for nf in range(n0, n1):

	print(f"processing file {nf-n0 + 1} of {n1-n0}")

	# Load configuration
	X, P = sbt.configuration_sorted(filenames[nf])

	# Compute mean square displacement
	dX[nf] = np.mean(np.sum((X - X0)**2, axis=1))

	# Compute orientation correlation
	for i in range(3):
		for j in range(3):
			dP[nf, i, j] = np.mean(P0[:, i] * P[:, j])

dX = dX[n0:n1]
dP = dP[n0:n1, :, :]

output = np.array([dX, dP[:, 0, 0], dP[:, 0, 1], dP[:, 0, 2], dP[:, 1, 0], dP[:, 1, 1], dP[:, 1, 2], dP[:, 2, 0], dP[:, 2, 1], dP[:, 2, 2]]) #reformat

# Save data
np.savetxt(f"{target}/meanSquareDisplacement.dat", output, delimiter=',')
