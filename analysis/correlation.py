
# Import modules
import numpy as np
import glob, sys, yaml
import sbt

from config import directories, target

# Get box size and volume fraction from command line input
try: nu, L = float(sys.argv[1]), int(sys.argv[2]) 
except: 
    print("Invalid input. Please provide volume fraction (nu) and box size (L).")
    sys.exit(1)

try: overwrite = int(sys.argv[3]) == 1
except: overwrite = False

# Locate filepath
source, target = sbt.locateData(nu, L, directories, f"{target}/correlation")

# Load in number of quadrature points
with open(f"{source}/CellConfig.yaml") as stream:
	numberOfQuadraturePoints = yaml.safe_load(stream)["cellHydro"][0]["numberOfQuadraturePoints"]

# Get list of files
allFiles = glob.glob(f"{source}/result/result*/SylinderDist_*.pvtp")

# Get number of files processed
processedFiles = glob.glob(f"{target}/corr*.dat") 
if overwrite: processedFiles = []

# Filter out processed files
processedNumbers = [int(file.split("-")[-1].split(".")[0]) for file in processedFiles]
allNumbers = [int(file.split("_")[-1].split(".")[0]) for file in allFiles]
filenames = [file for i, file in enumerate(allFiles) if allNumbers[i] not in processedNumbers]

# Number of files to process
Nfiles = len(filenames)
print(f"Found {Nfiles} unprocessed files.")

# Loop over all files
for counter, file in enumerate(filenames):

	# Get file number
	nf = int(file.split("_")[-1].split(".")[0])

	print(f"\rProcessing: {counter+1}/{Nfiles} (snapshot {nf})", end="", flush=True)

	# Load configuration
	X, P, Y, F = sbt.configuration(file, numberOfQuadraturePoints)

	# Compute correlation functions
	output = sbt.correlation(X, P, Y, F, L)

	# Save data
	np.savetxt(f"{target}/corr-{nf}.dat", output, delimiter=",")

print("\nAll done.")