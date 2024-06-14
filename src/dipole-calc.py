# get from Lucas, modified by Cong
import pandas as pd
import os, sys, re

bohr2ang = 0.529178

LiH_dipole_exp2 = 5.88

atomic_info_list = [
        [1.0079,1],
        [4.0026,2],
        [6.941, 3],
        [9.0122,4],
        [10.811,5],
        [12.011,6],
        [14.007,7],
        [15.999,8],
        [18.998,9],
]

atom_labels = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F']
atomic_info = pd.DataFrame(atomic_info_list, index=atom_labels, columns=['A', 'Z'])


def abs_pos(R, m1, m2, non_cm=False):
# f=3.5 for integrals_OH-0-2-bf59-noreorient
# reason for the formulae
# "prefactor of 3.5 was multiplied on the oxygen nuclear coordinate in Eq. (28)"
# J. Chem. Phys. 155, 154103 (2021)
# -> R1 = f R_parameter  m2/(m1+m2) (1) where R_parameter is a list in basis set generation
# R2 = - R_parameter m1/(m1+m2) (2)
# in computing dipole moment, 
# one starts from R in expec.t, that is R1 - R2
# R1 - R2 = R = R_parameter(fm2+m1)/(m1+m2) (3)
# from (3) we can solve R_parameter, thus obtain R_1 and R_2 as below
# namely, R_parameter = R * (m1 + m2) / (fm2+m1)
# R1 = f R_parameter  m2/(m1+m2) = f * R * (m1 + m2) / (fm2+m1) *  m2/(m1+m2) = R*f*((m2)/(m1 +f* m2))
# R2 = - R_parameter m1/(m1+m2) = -  R * (m1 + m2) / (fm2+m1) * m1/(m1+m2) = -R*((m1)/(f*m2 + m1))
#  f = 3.5
    f = 1.0 
    R1 = R*f*((m2)/(m1 +f* m2))
    R2 = -R*((m1)/(f*m2 + m1))

    if non_cm:
        R2 -= R1
        R1 -= R1
        print(R1,R2)

    return R1, R2

def debye2au(value):

    value = float(value)
    converted_value = value*0.393456

    return  float(converted_value)

def au2debye(value):

    value = float(value)
    converted_value = value/0.393456

    return float(converted_value)


def detect_expec_format(input_file):
    
    with open(input_file) as f_0:
        lines_0 = f_0.readlines()
        if 'nel' in lines_0[0]:
            flag_format = 1
        else:
            flag_format = 0
    
    return flag_format

def dipole_calc(outf=None):

    expect_file = 'expec.t'
    output_file = r'mcend.*.out'

    flag_format = detect_expec_format(expect_file)
    #print('flag_format',flag_format)

    Expec = pd.read_csv(expect_file, skiprows=flag_format, delimiter='\s+')
        
    R_expec = Expec.R.iloc[-1]
    z_expec = Expec.z.iloc[-1]

    for file_contents in os.listdir():
        outfile = re.match(output_file, file_contents)
        if outfile:
            outfile = outfile.group()
            break

    print('test',outfile)

    f = open(outfile, encoding="utf8", errors='ignore')
    flines = f.readlines()
    for line in flines:
        if 'Compound' in line:
            cmpd_line = line.split()
            cmpd = cmpd_line[1]
            sub_list = re.findall(r'[A-Z][a-z]*', cmpd)
            
        if 'electrons' in line:
            nel_line = line.split()
            #print('nel_line',nel_line)
            nel = nel_line[2]
            nel = int(nel)
            
            print('nel_line',nel_line, 'nel', nel)
            break

    print(sub_list)

    atm1 = sub_list[0]
    atm2 = sub_list[1]

    m1 = atomic_info.A[atm1]
    m2 = atomic_info.A[atm2]
    N1 = atomic_info.Z[atm1]
    N2 = atomic_info.Z[atm2]
    R1, R2 = abs_pos(R_expec, m1, m2)
    mu_e = z_expec  * -1 * nel
    mu_N = (N1*R1 + N2*R2)
    
    print(N1,N2,R_expec)
    print('m1,m2',m1,m2)
    print('R1,R2',R1,R2)
    print('mu_N',mu_N)
    print('mu_e',mu_e)
    
    mu = mu_N + mu_e

    return mu

dipole = au2debye(dipole_calc())

print('Dipole Moment: {: .3f} Debye'.format(dipole))


