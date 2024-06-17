#!/usr/bin/env python

#!/home/lucas/anaconda3/bin/python
#!/home/lea0105/miniconda3/envs/p4env/bin/python
#!/home/iulusoy/WORK/PROGRAMS/PSI4/miniconda/bin/python
#!/home/iulusoy/psi4conda/bin/python
# alter the path above
# should generalize at some point

import argparse
import psi4
from numpy import *
import os, sys, time, re
import pandas as pd
import writeints
from read_params import read_grid_parameters, read_poly_grid_parameters, write_molden_basis_string
from read_params import decomp_compound_name
from read_params import one_ints_write, two_eri_write
from read_params import read_coord
from gen_coord import *

from itertools import product
from pprint import pprint

#  bohr2Ang = 0.529178
bohr2Ang = psi4.constants.bohr2angstroms

parser = argparse.ArgumentParser(description='''
        Hi! I generate/orthogonalize the AO matrices for you using Psi4!
        ''')

epilog = '\nUsage: python do_psi4_ints_v2.py -f LiH -b cc-pvdz -r y -i y  '
parser.add_argument('-r', '--runtype', help='Moving basis. Leave out flag if basis is static (default=static)', required=False)
parser.add_argument('-b', '--basis', help='Which basis: cc-pvdz, 6-31G**, etc, quotes are not needed.', required=True)
parser.add_argument('-f', '--formula', help='what compound, to load rsp and grid_parameters', required=False)
parser.add_argument('-coord', '--sample_coordinate', help='sample internal coordinate', required=False)
parser.add_argument('-mc', '--charge', help='Molecular charge (default=0)', required=False)
parser.add_argument('-ms', '--multiplicity', help='Molecular multiplicity (default=1)', required=False)
parser.add_argument('-i', '--intsmethod', help='use python to write instead of fortran (default=no)', required=False)
parser.add_argument('-fc', '--fortcompiler', help='Fortran compiler to use (gfortran or ifort). Option deprecated', required=False)
# parser.add_argument('-gg', '--genguess',
#                     help='Generate an SCF values at the grid parameters. (Default=true)',
#                     required=False
#                     )

args = vars(parser.parse_args())

if args['runtype']:
    # unsupported
    runtype = 'moving'
else:
    runtype = 'static'

if args['basis']:
    basis = str(args['basis'])
else:
    raise Exception('Please specify a basis set')

if args['formula']:
    cmpd = str(args['formula'])
    input_coord_name = 'coord_' + cmpd 
    input_rsp_name   = 'rsp_' + cmpd
    input_grd_name   = 'grd_params_' + cmpd
    # generate atm list from the Z-matrix coord
    atm, type_list, type_list_total, link_list, value_list = read_coord(input_coord_name)
    
   # atm1, atm2 = decomp_compound_name(cmpd)
else:   
    raise Exception('Need to input a formula for molecule')

if args['charge']:
    mcharge = int(args['charge'])
else:
    mcharge = 0

if args['multiplicity']:
    m_mult = int(args['multiplicity'])
else:
    m_mult = 1


# gen_scf = bool(args['genguess'])
# don't really know why you'd need to do this without generating an initial guess
gen_scf = True


if args['intsmethod']:
    write_ints = 'python'
else:
    write_ints = 'fortran'
    
# atm from read_coord(input_coord_name) 
# input_coord_name: e.g., coord_LiH, prepared in advance, internal coordianate
print('before gen_poly_cartesian', atm, 'link_list', link_list, 'value_list', value_list)
# e.g., for HCN
#  ['C', 'H', 'N'] link_list [[], ['1'], ['1', '2']] value_list [[], ['1.064'], ['1.156', '180.']]
# from Z matrix to Cartisian, often used in generating ghost atoms
#gen_poly_cartesian(atm, link_list, value_list)
print('type_list', type_list, 'type_list_total', type_list_total)

maindb = pd.HDFStore('atomic_info.h5', 'r')
M = [None] * len(atm)
Z = [None] * len(atm) 
    
for i in range(len(atm)):
    M[i] = maindb.atomic.loc[atm[i], 'A']
    Z[i] = maindb.atomic.loc[atm[i], 'Z']
    
maindb.close()    
    
print('summary of M and Z', M, Z)  
print('input_rsp_name',input_rsp_name,'atm',atm, 'link_list',link_list,'value_list',value_list)
#sys.exit('test')

                                 #      rsp_LiH
                                 #                        ['Li', 'H']
                                 #                                     [[], ['1']]           
                                 #                                            [[], [4.0000082012776765]]        
cartesian_coordinates = gen_ghost_atom_positions(input_rsp_name, atm, link_list, value_list)    
print('cartesian_coordinates',cartesian_coordinates)

the_cutoff = 1e-8
threshd = 1.0e-1
dist_factor = 1.0

# load values from rsp, convert to angstroms
if runtype == 'static':
    #rgrid = loadtxt('rsp')*bohr2Ang
    # need to be generalized into multivariable cases
#    nrprimn, rmin, dr = read_grid_parameters(input_grd_name)
    nrprimn, rmin, dr = read_poly_grid_parameters(input_grd_name)
   # nrprimn as a list, e.g., [38], [38,39]

# or create spacing grid from parameters
elif runtype == 'moving':
    nrprimn, rmin, dr = read_grid_parameters()
    rmax = rmin + (nrprimn - 1)*dr
    #rgrid = linspace(rmin, rmax, nrprimn)*bohr2Ang




psi4.core.clean()
os.system('rm -f {}-psi4.out'.format(cmpd))
# all psi4 outputs will be written to this file
psi4.set_output_file("{}-psi4.out".format(cmpd), True)

# default orthog = symmetric

guess = 'SAD'
#guess = 'core'
scftyp = 'pk'
itermax = 255
#reference = 'rhf'
reference = 'uhf'
ps4_intoptions = {
    'basis':'{}'.format(basis),
    'print_mos':True,
    'print_basis':True,
    'guess':'{}'.format(guess),
    'reference':'{}'.format(reference),
    'scf_type':'{}'.format(scftyp),
    'maxiter':'{}'.format(itermax),
    'soscf':True,
    'SOSCF_START_CONVERGENCE': 1.0e-1      
#     's_orthogonalization':'canonical','symmetric',
#     's_tolerance':1e-8,
#     'e_convergence':1e-8,
#    'ints_tolerance':1e-12
}
psi4.set_options(ps4_intoptions)

t1 = time.time()

print('Generating Psi4 integrals for {} grid points'.format(len(cartesian_coordinates)))


for k in arange(1, len(cartesian_coordinates)+1):
    # if this is the first run, write all integrals out
    write_all = bool(k == 1)

    time_init = time.time()
    geom_str = ''
    geom_str += 'nocom \n'   
    geom_str += 'noreorient \n'    
    #> write coordinates so ghost atoms are on 3 locations and real atoms on 1
    #> then iterate through every combination
    #> added a fix to get generate standard dipole, was backwards before due to
    #> a weird oversight
#     if runtype == 'static':
    for i in range(len(cartesian_coordinates)):
        if i + 1 == k:
            # new way
#             geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm2,  posM2[i])
#             geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm1, -posM1[i])
            # old way
#             geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm1, -posM1[i])
#             geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm2, -posM2[i])
#           geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm1, posM1[i])
#           geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm2, posM2[i])
            for m in range(len(atm)):
              geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm[m], cartesian_coordinates[i][m][2])            
            
        else:
            # new way
#             geom_str += '@{}   0.0  0.0  {: 9.7f}\n'.format(atm2,  posM2[i])
#             geom_str += '@{}   0.0  0.0  {: 9.7f}\n'.format(atm1, -posM1[i])
            # old way
#            geom_str += '@{}   0.0  0.0  {: 9.7f}\n'.format(atm1, posM1[i])
#            geom_str += '@{}   0.0  0.0  {: 9.7f}\n'.format(atm2, posM2[i])
            for m in range(len(atm)):
              geom_str += '@{}   0.0  0.0  {: 9.7f}\n'.format(atm[m], cartesian_coordinates[i][m][2])   
              
#     elif runtype == 'moving':
#         i = k - 1
#         geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm1,  posM1[i])
#         geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm2, -posM2[i])

    geom_str += 'symmetry c1\n'
    print('geom_str', geom_str, type(geom_str))
    

    mol = psi4.geometry(geom_str)
    print('mol', mol)
    #print('psi',  mol.print_out()  )
    #> set the charge
    mol.set_molecular_charge(mcharge)

   # print('geom',geom_str)

    print('m_mult',m_mult)

    mol.set_multiplicity(m_mult) 

    #> do scf calculation and call integrals
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS'))
    mints = psi4.core.MintsHelper(wfn.basisset())
    
    Etot, wfn = psi4.energy('SCF', return_wfn=True)



    print('Etot',Etot)

    os.system('rm -f {}_r{}.molden'.format(cmpd, k))
    psi4.molden(wfn, '{}_r{}.molden'.format(cmpd, k))

    T = asarray(mints.ao_kinetic(), order='F')
    nbf = int(shape(T)[0])
    V = asarray(mints.ao_potential(), order='F')
    S = asarray(mints.ao_overlap(), order='F')
    r = mints.ao_dipole()

    print('write_all', write_all)


    if write_all:
        Vee = asarray(mints.ao_eri(), order='F').reshape(nbf, nbf, nbf, nbf)

    print('Writing all integrals')
    H = T + V
    x = array(r[0], order='F')
    y = array(r[1], order='F')
    z = array(r[2], order='F')
    #> write in tridiagonal form
    il1 = tril_indices(nbf)

    print('len',len(T),len(H))
    print('shape', shape(T),shape(H) )
    print('T',T[0,0],T[0].shape)
    print('H',H[0],H[0].shape)


    H = H[il1]
    S = S[il1]
    T = T[il1]
    x = x[il1]
    y = y[il1]
    z = z[il1]
    nrelem = len(H)
    print('after triindices', T.shape,T[0])
    #breakpoint()
    print('writeints', write_ints)
   
    #> if the f2py generation isn't successful use back up way with python
    #> substantially slower but with the 8-fold symmetry it isn't as bad
    #> still best to avoid for larger systems
    if write_ints == 'python':
        intfile = open('fort.999_{:1d}'.format(k), 'w')


        print('intfile',intfile)

        hts_mat = one_ints_write(H, T, S, nbf)
        intfile.write(hts_mat)
        xyz_mat = one_ints_write(x, y, z, nbf)
        intfile.write(xyz_mat)
        if write_all:
            two_ints, nrteint = two_eri_write(Vee, nbf)
            intfile.write(two_ints)
        intfile.close()

    #> Call fortran program to write files, because python is way too slow at writing to output
    elif write_ints == 'fortran':

        print('wrtie_all', write_all)

        #> write out the two-electron integrals with the 8-fold symmetry
        if write_all:

            print('rij',write_all)
            print('output file') 
            print('fort.999_{:1d}'.format(k))

           # intfile = open('fort.999_{:1d}'.format(k), 'w')
           # print('intfile',intfile)



            n = nbf
            rijindex = []
            rij = []
            for ix in range(n):
                for jx in range(ix+1):
                    ij = (ix*(ix + 1)//2 + jx)
                    for kx in range(n):
                        for lx in range(kx+1):
                            kl = (kx*(kx + 1)//2 + jx)
                            if ij >= kl:
                                if abs(Vee[ix, jx, kx, lx]) > the_cutoff:
                                    rijindex.append([ix+1, jx+1, kx+1, lx+1])
                                    rij.append(Vee[ix, jx, kx, lx])
            rijindex = array(rijindex, order='F', dtype=int)
            rijindex = rijindex.T
            rij = array(rij, order='F')
            nrteint = len(rij)

            print('integral writer file', 'fort.999_{:1d}'.format(k) ) 
#          subroutine integral_writer(H, T, S, x, y, z, rij, rijindex, outfilename,  nrelem, nrtwoeri)
            writeints.integral_writer(H, T, S, x, y, z, rij, rijindex, 'fort.999_{:1d}'.format(k), nrelem, nrteint)
        else:
            print('else',write_all)
            writeints.hm_writer(H, T, S, x, y, z, 'fort.999_{:1d}'.format(k), nrelem)

    #> only need to do this once
    print('write_all', write_all)

    if write_all:
        print('writing to data.inp')
        dfile = open('data.inp', 'w')
#         dfile.write('{}\n'.format(nbf))
#         dfile.write('{}\n'.format(int(nrteint)))
#         dfile.write('{}\n'.format(int(len(cartesian_coordinates))))

        dfile.write('{}\n{}\n{}\n'.format(nbf, int(nrteint), int(len(cartesian_coordinates))))

        dfile.close()

    if write_all:
        print('Main Pass completed in {:7.2f} s'.format(time.time() - time_init))
    else:
        print('Pass {:3d} of {:3d} completed'.format(k, len(cartesian_coordinates)))

if runtype == 'moving':
    scfguess.close()

t2 = time.time()
write_time = t2 - t1
print('All integrals written to file. Total Runtime: {:7.2f} s'.format(write_time))

#> finder script
mbasis_name = write_molden_basis_string(cmpd, basis)
print('SCF calculation output written to {}-psi4.out'.format(cmpd))




def gen_scf_guess(atm, M, cmpd, link_list, value_list, basis, mcharge, dist_factor, psiopts=ps4_intoptions):

#     nrprimn, rmin, dr = read_grid_parameters()
    dist_factor_2 = 0.75
    #dist_factor_2 = 1.00
    
    n_dof = len(nrprimn)
    print('n_dof',n_dof)
    
    for i in range(n_dof):
    
      rmax = rmin[i] + (nrprimn[i] - 1)*dist_factor_2*dr[i]
   # rgrid = linspace(rmin, rmax, nrprimn)*bohr2Ang
      rgrid = linspace(rmin[i], rmax, nrprimn[i])#*bohr2Ang
   
      print('dist_factor', dist_factor)
      print('rgrid', rgrid)
    # need to be revised such that one series one SPF
    # assume dof order
      gen_grd_file('grd_params_' + cmpd, 'grd_' + cmpd)
      #input: a.u.
      #output: angstrom
      #cartesian_coordinates = gen_nrgrid_positions('grd_' + cmpd, atm, link_list, value_list)
      cartesian_coordinates = gen_nrgrid_positions_v2(atm, link_list, value_list, rgrid, i)
     # cartesian_coordinates = gen_nrgrid_positions_v2(rgrid, i)
      
      print('cartesian_coordinates!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', len(cartesian_coordinates))
    
    
    
  #  sys.exit()

# NEED REGENERATE!!    
#   posM1 = rgrid*dist_factor*M2/(M1 + M2)
#   posM2 = -rgrid*M1/(M2 + M1)
#   posM2 = rgrid*M1/(M2 + M1)

      psi4.set_options(psiopts)
      scfguess = open('scf-guess-{}.dat'.format(cmpd), 'w')
    

   # for k in arange(1, len(rgrid)+1):
      for k in arange(1, len(cartesian_coordinates)+1):
          time_init = time.time()
          geom_str = ''

          i = k - 1
        #geom_str += ' {}   0.0    0.0   {: 9.7f}\n'.format(atm2, posM2[i])
        #geom_str += ' {}   0.0    0.0   {: 9.7f}\n'.format(atm1, -posM1[i])
        
          for m in range(len(atm)):
              geom_str += ' {}   0.0  0.0  {: 9.7f}\n'.format(atm[m], cartesian_coordinates[i][m][2])         
        
        
          geom_str += 'symmetry c1\n'


          print('geom_str',geom_str)

          mol = psi4.geometry(geom_str)

          mol.set_molecular_charge(mcharge)
          wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS'))
          mints = psi4.core.MintsHelper(wfn.basisset())

          Etot = psi4.energy('SCF')
          print('Etot',Etot)
        # could use full-ci for energy guess but doesn't change all that much
#         Etot = psi4.energy('FCI')

          scfguess.write('{: 23.16e}\n'.format(Etot))

    scfguess.close()
    

if gen_scf:
    print('Generating SCF initial guess for {} at each grid point'.format(cmpd))
    gen_scf_guess(atm, M, cmpd, link_list, value_list, basis, mcharge, dist_factor)
    print('Done.')

os.system('rm -f *.mod')

    

# if not args['fortcompiler'] or args['fortcompiler'] == 'gfort' or args['fortcompiler'] == 'gfortran':
#     os.system("gfortran write2string.f90 generate_int_v2.0.f90 -llapack -fopenmp -lblas -O3 -o generate_int_v2.0.o")

#> should be part of make script, whenever I figure out how to add all that together properly
# elif args['fortcompiler'] == 'ifort':
# #     os.system("ifort write2string.f90 generate_int_v2.0.f90 -mkl -qopenmp -O3 -fp-model precise -o generate_int_v2.0.o")
# #     os.system("ifort write2string.f90 generate_int_v2.0.f90 -mkl -qopenmp -O3 -o generate_int_v2.0.o")
#     os.system("export FC=ifort")
#     os.system("make generate")

# elif runtype == 'moving':
#     rspn, rmin, rmax = read_grid_parameters()
#     print('Using moving basis with {} grid points'.format(rspn))
#     raise Exception('Option not yet supported')


f = open('threshd.dat', 'w')
f.write('{: 23.16e}'.format(threshd))
f.close()
os.system('./generate_int_v2.0.o')
# os.system('rm -f threshd.dat')


# get number of non-zero two-electron integrals and num of basis functions kept
f = open('basis-purge.dat', 'r')
for lines in f.readlines():
    if re.match(r'^\s*\d+ Nonzero two-electron', lines):
        nrtwoeint = int(re.findall(r'\d+', lines)[0])
        print('nint2e = ', nrtwoeint)
#         break
    if re.match(r'^\s*\d+ basis functions remain', lines):
        nrprime = int(re.findall(r'\d+', lines)[0])
        print('basis functions = ', nrprime)
#         break
f.close()



rspn = len(cartesian_coordinates)
hmat_file = '{}-hmat-r{:02d}-{:02d}.dat'.format(cmpd, 1, rspn)
ints_file = '{}-ints-{:02d}.dat'.format(cmpd, 1)
scf_file = 'scf-guess-{}.dat'.format(cmpd)

nr2file = open('{}'.format(basis), 'w')
nr2file.write('nint2e = {}'.format(nrtwoeint))
nr2file.close()

info_file = open('info-{}.dat'.format(cmpd), 'w')
info_file.write('nint2e = {}\n'.format(nrtwoeint))
info_file.write('nrprime = {}\n'.format(nrprime))
info_file.write('nrprimn = {}\n'.format(nrprimn))
info_file.write('nrensp = {}\n'.format(rspn))
info_file.write('rmin = {}\n'.format(rmin))
info_file.write('dr = {}'.format(dr))
info_file.close()


os.rename('data_hmat.asc', hmat_file)
os.rename('data_out.asc', ints_file)

ints_dir = 'integrals_{}-bf{}'.format(cmpd, nrprime)
if not os.path.isdir(ints_dir):
    os.mkdir(ints_dir)

os.system('cp smateigvals.n {}'.format(ints_dir))
os.system('mv smat_xmat.dat {}'.format(ints_dir))
# need to be revised to rsp_ -> rsp
#os.system('cp rsp {}/rsp_{}'.format(ints_dir, cmpd))
os.system('cp rsp_{} {}/rsp_{}'.format(cmpd, ints_dir, cmpd))
os.system('mv {} {}'.format(hmat_file, ints_dir))
os.system('mv {} {}'.format(ints_file, ints_dir))
os.system('mv {} {}'.format(basis, ints_dir))
os.system('mv {} {}'.format('info-{}.dat'.format(cmpd), ints_dir))
os.system('mv {} {}'.format(mbasis_name, ints_dir))
if gen_scf:
    os.system('mv {} {}'.format(scf_file, ints_dir))

# move non-vital files to additional folder to reduce clutter
aux_dir = '{}_aux_basisfiles'.format(cmpd)
if not os.path.isdir(aux_dir):
    os.mkdir(aux_dir)

os.system('mv fort.999_* {}'.format(aux_dir))
os.system('mv {}_r*.molden {}'.format(cmpd, aux_dir))
os.system('mv smateigvals.n {}'.format(aux_dir))
os.system('mv one {}'.format(aux_dir))
os.system('mv {}-psi4.out {}'.format(cmpd, aux_dir))



