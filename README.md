# CNMFw3

CNMFw3 (CNMF for wide-field 3D image) is a Constrained NMF-based neural activity extraction algorithm for 3-dimiensional calcium imaging movie.

## Usage

Load movie data as 4d array named Y (XYZT) and run runCNMFw3.m


## Dependencies

CNMFw3 has the following dependencies:
- MATLAB (2018a or later)
- CNMF (https://github.com/flatironinstitute/CaImAn-MATLAB)
- ordfilt3 (https://www.mathworks.com/matlabcentral/fileexchange/22044-ordfilt3) 

## Installation

1. Download or clone this repository and add its path to your MATLAB path.

2. Download and install the dependent softwares listed above.

## Example

demoCNMFw3: this demo runs CNMFw3 with simulated data (100 x 100 x 10 voxels, 1000 frames).

>>[filters,traces]=demoCNMFw3;

This demo took <100 seconds on a desktop PC with Intel Core i7-4790 and 32GB RAM.

## License

This project is covered under the MIT license.
