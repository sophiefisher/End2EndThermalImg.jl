import argparse
import grcwa as rc
import numpy as np

def get_transmission(freq, pillar_width, pillar_height, pillar_epsilon, unit_cell_length, substrate_epsilon, nG):
    # nG is number of fourier components in RCWA (truncation order)
    L1 = [unit_cell_length,0]
    L2 = [0,unit_cell_length]
    N_grid_points = 10000 # number of grid points for the patterned layers
    obj = rc.obj(nG, L1, L2, freq, 0.0, 0.0, verbose = 0 )
    obj.Add_LayerUniform(3.0, substrate_epsilon)
    obj.Add_LayerGrid(pillar_height, N_grid_points, N_grid_points)
    obj.Add_LayerUniform(3.0, 1.0)  
    obj.Init_Setup()
    epgrid = np.zeros((N_grid_points, N_grid_points))
    boundary = unit_cell_length / 2.0
    xarr = np.linspace(-boundary, boundary, N_grid_points)
    yarr = np.linspace(-boundary, boundary, N_grid_points)
    x, y = np.meshgrid(xarr, yarr, indexing='ij' )
    ind = np.logical_and(abs(x) < pillar_width / 2.0, abs(y) < pillar_width / 2.0)
    epgrid[ind] = 1.0
    epgrid = (pillar_epsilon - 1.0) * epgrid + 1.0
    obj.GridLayer_geteps(epgrid.flatten())
    obj.a0 = np.zeros(2*obj.nG, dtype=complex)
    obj.a0[0] = 1.0
    obj.bN = np.zeros(2*obj.nG, dtype=complex)
    aN, b0 = obj.GetAmplitudes(2, 0.1)
    transmission = aN[0]/obj.a0[0] * np.sqrt(np.sqrt(substrate_epsilon))
    return transmission

def main():
    parser = argparse.ArgumentParser(description="Compute transmission using RCWA.")
    parser.add_argument("freq", type=float, help="Frequency of incident wave")
    parser.add_argument("pillar_width", type=float, help="Width of the pillar")
    parser.add_argument("pillar_height", type=float, help="Height of the pillar")
    parser.add_argument("pillar_epsilon", type=float, help="Dielectric constant of the pillar")
    parser.add_argument("unit_cell_length", type=float, help="Unit cell length")
    parser.add_argument("substrate_epsilon", type=float, help="Dielectric constant of the substrate")
    parser.add_argument("nG", type=int, help="Number of Fourier components")
    
    args = parser.parse_args()
    
    transmission = get_transmission(
        args.freq, args.pillar_width, args.pillar_height,
        args.pillar_epsilon, args.unit_cell_length, args.substrate_epsilon, args.nG
    )
    
    print(f"{transmission.real} + {transmission.imag}im")

if __name__ == "__main__":
    main()