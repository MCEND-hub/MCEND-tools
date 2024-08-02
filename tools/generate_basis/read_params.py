import sys
import re

# auxilary functions for the integral routines

# params file structure
#  spacing and order do not matter
#      nrprimn = 128
#      rmin = 0.68030135
#      dr = 0.01889726


def read_grid_parameters():
    parf = open("grd_params", "r")
    flt_regex = "[-+]?([0-9]*\.[0-9]+|[0-9]+)"
    nrprimn = False
    rmin = False
    dr = False
    for line in parf.readlines():
        if "nrprimn" in line:
            nrprimn = re.findall(flt_regex, line)
            nrprimn = int(nrprimn[0])
        if "rmin" in line:
            rmin = re.findall(flt_regex, line)
            rmin = float(rmin[0])
        if "dr" in line:
            dr = re.findall(flt_regex, line)
            dr = float(dr[0])
    # if any of the three parameters are still set as false
    # then something didn't go right
    if not nrprimn or not rmin or not dr:
        print("params file was not read correctly please check it")
        sys.exit()
    else:
        pass

    parf.close()

    return nrprimn, rmin, dr


# I know it is not a good practice
# separate then combine
# N.B. need to be merged
def read_poly_grid_parameters(input_file):
    parf = open(input_file, "r")
    flt_regex = "[-+]?([0-9]*\.[0-9]+|[0-9]+)"
    nrprimn = False
    rmin = False
    dr = False
    for line in parf.readlines():
        if "nrprimn" in line:
            nrprimn = re.findall(flt_regex, line)
            nrprimn = [int(x) for x in nrprimn]
        # nrprimn = int(nrprimn[0])

        if "rmin" in line:
            rmin = re.findall(flt_regex, line)
            rmin = [float(x) for x in rmin]
        #       rmin = float(rmin[0])
        if "dr" in line:
            dr = re.findall(flt_regex, line)
            dr = [float(x) for x in dr]

    #        dr = float(dr[0])
    # if any of the three parameters are still set as false
    # then something didn't go right
    if not nrprimn or not rmin or not dr:
        print("params file was not read correctly please check it")
        sys.exit()
    else:
        pass

    parf.close()

    return nrprimn, rmin, dr


# read internal coordinate
# since it's bit difficult to specify quantum dof from cartesian coordinate
# e.g. https://github.com/psi4/psi4/blob/master/samples/opt6/input.dat
#   O
#   H 1 1.0
#   H 1 1.0 2 160.0
# input: Li
# H 1 1.596
# output: [[],[1.596]]
# HCN
# [[], ['1'], ['1', '2']]
# [[], ['1.064'], ['1.156', '180.']]
def read_coord(input_file):

    atm = []
    # atom list
    type_list = []
    # dist, angle, dihedral
    type_list_total = []
    # [[1],[2],[3]]
    link_list = []
    # [to 2,3,4...]
    value_list = []
    # [[1],[2],[3]]
    # not sure if this is the best approach
    #    mol_coord = open('coord','r')
    mol_coord = open(input_file, "r")
    for n, line in enumerate(mol_coord.readlines()):

        print(n, line)
        temp_type_list = []
        temp_link_list = []
        temp_value_list = []

        if len(line.strip()) == 0:
            continue

        if len(line) > 1:
            elements = line.split()
        else:
            elements = line

        print("elements", elements)
        atm.append(elements[0])

        for i, element in enumerate(elements):
            if i == 0:
                continue

            if i in [1, 3, 5]:
                temp_link_list.append(element)

            if i in [2]:
                temp_type_list.append("dist")
                temp_value_list.append(element)

            if i in [4]:
                temp_type_list.append("angle")
                temp_value_list.append(element)

            if i in [6]:
                temp_type_list.append("dihedral")
                temp_value_list.append(element)

        type_list.append(temp_type_list)
        link_list.append(temp_link_list)
        value_list.append(temp_value_list)

    for type_elements in type_list:
        for type_element in type_elements:
            type_list_total.append(type_element)

    print(atm)
    print(type_list)
    print(type_list_total)
    print(link_list)
    print(value_list)

    return atm, type_list, type_list_total, link_list, value_list


def one_ints_write(A, B, C, nbf):
    int_str = ""
    for i in range(len(A)):
        int_str += "{: 24.16e} {: 24.16e} {: 24.16e}\n".format(A[i], B[i], C[i])
    return int_str


def two_eri_write(eri, nbf):
    n = nbf
    int_str = ""
    nrtwoeint = 0

    for ix in range(n):
        for jx in range(ix + 1):
            ij = ix * (ix + 1) // 2 + jx
            for kx in range(n):

                # print('kx',kx)

                for lx in range(kx + 1):

                    # print('lx',lx)

                    kl = kx * (kx + 1) // 2 + jx
                    if ij >= kl:
                        if abs(eri[ix, jx, kx, lx]) > 1e-15:
                            int_str += "{:>5} {:>5} {:>5} {:>5} {: 24.16e}\n".format(
                                ix + 1, jx + 1, kx + 1, lx + 1, eri[ix, jx, kx, lx]
                            )
                            nrtwoeint += 1
    #                     ij = int(str(ix) + str(jx))
    #                     kl = int(str(kx) + str(lx))
    return int_str, nrtwoeint


def decomp_compound_name(usr_inp):
    element = re.findall(r"[A-Z][a-z]*|\d+|\(|\)", usr_inp)
    if len(element) != 2:
        print("Please enter a diatomic molecule")
        sys.exit()

    else:
        atm1 = element[0]
        if element[1].isalpha():
            atm2 = element[1]
        elif element[1].isdigit:
            atm2 = element[0]
        return atm1, atm2


def write_molden_basis_string(cmpd, basisset):
    #     f = open('{}_r2.molden'.format(cmpd), 'r')
    f = open("{}_r1.molden".format(cmpd), "r")
    mlines = f.readlines()
    mbasis = ""

    for i, line in enumerate(mlines):
        if "[GTO]" in line:
            for ix in range(i, len(mlines)):
                if "[MO]" not in mlines[ix]:
                    mbasis += mlines[ix]
                else:
                    break
    f.close()

    mfile_name = "{}-{}.mbasis".format(cmpd, basisset)
    mfile = open(mfile_name, "w")
    print(mbasis, file=mfile)
    mfile.close()
    return mfile_name


def print_matrix(mat, header, print_file=None):

    cols = mat.shape[0]
    rows = mat.shape[1]
    str_mat = "{:^{width}}\n".format(header, width=rows * 10)
    for i in range(cols):
        str_mat += "{}\n".format(" ".join("{:> 12.6f}".format(x) for x in mat[i, :]))
    print(str_mat, file=print_file)
