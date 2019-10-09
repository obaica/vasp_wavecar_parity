# vasp_WAVECAR_parity

## Author

 Wang Dinghui
 
 DZ1722031@smail.nju.edu.cn
 
 Nanjing University

## Usage:

It already has a precompiled excutable program "parity.x", or you may compile by

> ifort -o parity.x wavetrans.f90 parity.f90 main.f90 -assume byterecl

and running by
> /dir_to_program/parity.x

### Please make sure the following points:
- Setting LWAVE=.TRUE. in your INCAR and remaining WAVECAR.
- Target k point included in your KPIONTS
- Your system DOSE have Inversion Symmetry.
- Exact Inverstion Symmetry Coordinate corresponding to POSCAR.

### tips
- You may run a non-scf calculation which has few kpoints ( target K included of course!) after scf calculation.
- About INVERSION SYMMETRY, usually it's not locate in (0,0,0). So you should carefully check it's coordinate according to POSCAR or you'll get wrong result.
- For the case considering soc, if the system has time reversal symmetry then each band is doubly degenerate and only one need to take into count. 

### warning:
- Now this code is not tested in magnetic system yet!

## Reference
WAVETRANS.F90 code modified from https://www.andrew.cmu.edu/user/feenstra/wavetrans/  by R. M. Feenstra and M. Widom Department of Physics, Carnegie Mellon University, Pittsburgh, PA 15213
