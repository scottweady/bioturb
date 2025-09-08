
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
source, target = sbt.locateData(nu, L, directories, f"{target}/correlation")

# Load in number of quadrature points
with open(f"{source}/CellConfig.yaml") as stream:
	numberOfQuadraturePoints = yaml.safe_load(stream)["cellHydro"][0]["numberOfQuadraturePoints"]

# Get number of files processed
processedFiles = glob.glob(f"{target}/corr*.dat") 
Nprocessed = len(processedFiles)

# Get list of files and sort
filenames = glob.glob(f"{source}/result/result*/SylinderDist_*.pvtp")
filenames = natsorted(filenames)
Nfiles = len(filenames)

# Loop over all files
for nf in range(Nprocessed, Nfiles):

	print(f"Processing file {nf+1} of {Nfiles}")

	# Load configuration
	X, P, Y, F = sbt.configuration(filenames[nf], numberOfQuadraturePoints)

	# Compute correlation functions
	output = sbt.correlation(X, P, Y, F, L)

	# Save data
	np.savetxt(f"{target}/corr-{nf}.dat", output, delimiter=",")
