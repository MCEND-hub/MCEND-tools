import psi4
import numpy as np
import copy

# import pprint
from itertools import product

bohr2Ang = psi4.constants.bohr2angstroms


def gen_ghost_atom_positions(
    input_rsp_file_name, atm, link_list, value_list
):  # ,type_list_total):
    #                              rsp_HCN
    # e.g.,  ['C', 'H', 'N'] [[], ['1'], ['1', '2']] [[], ['1.064'], ['1.156', '180.']]
    # a.u. input
    print(
        "gen_ghost_atom_positions input",
        input_rsp_file_name,
        atm,
        link_list,
        value_list,
    )

    with open(input_rsp_file_name) as rsp_content:
        rsp_lines = rsp_content.readlines()
        pos_list = []
        for rsp_line in rsp_lines:
            if rsp_line[0] == "#":
                continue

            if rsp_line in ["\n", "\r\n"]:
                continue

            # temp_pos_list = []
            rsp_list = rsp_line.split()[0:-3]
            rsp_list = [float(i) for i in rsp_list]
            rsp_list = np.array(rsp_list)
            print("rsp_list", rsp_list)

            rgrid = bohr2Ang * rsp_list
            rgrid = list(rgrid)
            print("rgrid in gen_ghost_atom_positions", rgrid)

            connection_list = rsp_line.split()[-3:-1]
            connection_list = [int(i) for i in connection_list]
            print("connection_list", connection_list)

            type_list = rsp_line.split()[-1:]
            print("type_list", type_list)

            if type_list[0] == "dist":
                pos_list.append(rgrid)
            else:
                print("not supported yet")

        # list of rgrids, intented to include angular variables, so far not
        # [ [...], [...] ]
        print("pos_list", pos_list)

        # generate combination of grids, i.e., atomic positions
        # sort of flat
        list_list_gen = []
        for i in range(len(pos_list)):
            list_list_gen_temp = []
            for j in range(len(pos_list[i])):
                list_list_gen_temp.append(j)

            list_list_gen.append(list_list_gen_temp)

        # print('list_list_gen', list_list_gen)
        list_list = list(product(*list_list_gen))
        # list_list [(0,), (1,), (2,), (3,)]
        print("list_list", list_list)

        output_rsp_list = []
        for l in list_list:
            temp_output_rsp_list = []
            # print('l',l)
            # print(l[0],l[1],type(l[0]),type(l[1]))
            for n, index in enumerate(l):
                # print( 'pos', pos_list[n][index] )
                # weird index
                # [[0.8000016402555353, 1.6000032805110707, 2.8000057408943735, 4.0000082012776765]]
                temp_output_rsp_list.append(pos_list[n][index])

            output_rsp_list.append(temp_output_rsp_list)
        #          else:
        #              print('not supported')
        # e.g., output_rsp_list 4 [[0.8000016402555353], [1.6000032805110707], [2.8000057408943735], [4.0000082012776765]]
        print("output_rsp_list", len(output_rsp_list), output_rsp_list)

        gen_cartesian_input_list = []

        # gathering z-matrix
        for i, output_rsp in enumerate(output_rsp_list):
            # print('i',i)
            temp_list = []
            temp_list.append(atm)
            temp_list.append(link_list)

            temp_value_list = copy.deepcopy(value_list)
            # e.g., temp_value_list [[], ['1.596']]
            print("temp_value_list", temp_value_list)
            print("output_rsp", output_rsp)

            for j, value in enumerate(temp_value_list):
                if j == 0:
                    continue
                else:
                    # for (k,rsp_entry) in enumerate(output_rsp):
                    # provide temp_value from rsp lists
                    temp_value_list[j][0] = output_rsp[j - 1]

            temp_list.append(temp_value_list)
            # print('temp_list', temp_list)
            print(
                "atm", atm, "link_list", link_list, "temp_value_list", temp_value_list
            )
            temp_cartesian = gen_poly_cartesian(atm, link_list, temp_value_list)

            # gen_z_matrix_input_list.append(temp_list)
            #  gen_cartesian_input_list.append(temp_list)
            gen_cartesian_input_list.append(temp_cartesian)
    # gen_cartesian_input_list = adjust_coordinates(gen_cartesian_input_list)
    gen_cartesian_input_list = adjust_coordinates(gen_cartesian_input_list)
    # gen_cartesian_input_list = np.array(gen_cartesian_input_list)
    print(
        "gen_cartesian_input_list",
        gen_cartesian_input_list,
        type(gen_cartesian_input_list),
    )
    # generate z-matrix
    return gen_cartesian_input_list

    # sys.exit("test")


# def gen_grd_file(input_file,output_file):
def gen_grd_file(input_file, output_file):
    print("gen_grd_file")
    with open(input_file) as grd_content:
        grd_lines = grd_content.readlines()

        for n, grd_line in enumerate(grd_lines):
            grd_line = grd_line.strip("\n")
            elements = grd_line.split(" ")
            print(elements)
            elements = list(filter(None, elements))
            print(elements)

            if elements[0] == "nrprimn":
                nrprimn_list = elements[2:]
                nrprimn_list = [int(i) for i in nrprimn_list]
            if elements[0] == "rmin":
                rmin_list = elements[2:]
                rmin_list = [float(i) for i in rmin_list]
            if elements[0] == "dr":
                dr_list = elements[2:]
                dr_list = [float(i) for i in dr_list]
            if elements[0] == "link":
                link_list = elements[2:]
                link_list = [int(i) for i in link_list]
            if elements[0] == "type":
                type_list = elements[2:]
                # type_list = [int(i) for i in dr_list]

        print(nrprimn_list, rmin_list, dr_list, link_list, type_list)
        output_grd_list = []

        for n in range(len(nrprimn_list)):
            # for i in range(nrprimn
            rmin_temp = rmin_list[n]
            nrprimn_temp = nrprimn_list[n]
            dr_temp = dr_list[n]
            rmax_temp = rmin_temp + (nrprimn_temp - 1) * dr_temp
            #            rgrid_temp = np.linspace(rmin_temp, rmax_temp, nrprimn_temp)*bohr2Ang
            rgrid_temp = np.linspace(rmin_temp, rmax_temp, nrprimn_temp)  # *bohr2Ang

            print("rgrid_temp", type(rgrid_temp), rgrid_temp)

            rgrid_temp = list(rgrid_temp)

            for link in link_list:
                rgrid_temp.append(link)

            for t in type_list:
                rgrid_temp.append(t)

            # print('rgrid_temp', type(rgrid_temp), rgrid_temp)
            rgrid_temp_str = " ".join(str(v) for v in rgrid_temp)
            # print('rgrid_temp_str', type(rgrid_temp_str), rgrid_temp_str)

            output_grd_list.append(rgrid_temp_str)

        # print('output_grd_list',output_grd_list)

        output_str = "\n".join(output_grd_list)
        print("output_grd_list", output_str)
        print("output_file", output_file)
        f = open(output_file, "w")
        f.write(output_str)
        f.close()


def gen_nrgrid_positions(
    input_grd_file_name, atm, link_list, value_list
):  # ,type_list_total):
    # e.g.,  ['C', 'H', 'N'] [[], ['1'], ['1', '2']] [[], ['1.064'], ['1.156', '180.']]
    print(
        "gen_nrgrid_positions input",
        "atm",
        atm,
        "link_list",
        link_list,
        "value_list",
        value_list,
    )

    with open(input_grd_file_name) as grd_content:
        grd_lines = grd_content.readlines()
        pos_list = []
        for grd_line in grd_lines:
            if grd_line[0] == "#":
                continue

            if grd_line in ["\n", "\r\n"]:
                continue

            print("grd_line", grd_line)

            grd_list = grd_line.split()[0:-3]
            grd_list = [float(i) for i in grd_list]
            grd_list = np.array(grd_list)
            print("grd_list", grd_list)

            rgrid = bohr2Ang * grd_list
            rgrid = list(rgrid)
            print("rgrid in gen_ghost_atom_positions", rgrid)

            connection_list = grd_line.split()[-3:-1]
            connection_list = [int(i) for i in connection_list]
            print("connection_list", connection_list)

            type_list = grd_line.split()[-1:]
            print("type_list", type_list)

            if type_list[0] == "dist":
                pos_list.append(rgrid)
            else:
                print("not supported yet")

        # list of rgrids, intented to include angular variables, so far not
        print("pos_list", pos_list)
        # e.g.,  [[0.49999996680528824, 0.5999753534260497, 0.6999507400468112

        # generate combination of grids, i.e., atomic positions
        list_list_gen = []
        for i in range(len(pos_list)):
            list_list_gen_temp = []
            for j in range(len(pos_list[i])):
                list_list_gen_temp.append(j)

            list_list_gen.append(list_list_gen_temp)

        # print('list_list_gen', list_list_gen)
        # no need for combination
        # one series one SPF

        list_list = list(product(*list_list_gen))
        print("list_list", list_list)

        output_rsp_list = []
        for l in list_list:
            temp_output_rsp_list = []
            for n, index in enumerate(l):
                temp_output_rsp_list.append(pos_list[n][index])

            output_rsp_list.append(temp_output_rsp_list)

        print("output_rsp_list", len(output_rsp_list), output_rsp_list)

        gen_cartesian_input_list = []

        # gathering z-matrix
        for i, output_rsp in enumerate(output_rsp_list):
            # print('i',i)
            temp_list = []
            temp_list.append(atm)
            temp_list.append(link_list)

            temp_value_list = copy.deepcopy(value_list)

            # print('temp_value_list',temp_value_list)
            # print('output_rsp', output_rsp)

            for j, _ in enumerate(temp_value_list):
                if j == 0:
                    continue
                else:
                    # for (k,rsp_entry) in enumerate(output_rsp):
                    # provide temp_value from rsp lists
                    temp_value_list[j][0] = output_rsp[j - 1]

            temp_list.append(temp_value_list)
            # print('temp_list', temp_list)
            # print('atm',atm, 'link_list', link_list, 'temp_value_list', temp_value_list)
            print("temp_value_list", temp_value_list)
            temp_cartesian = gen_poly_cartesian(atm, link_list, temp_value_list)

            # gen_z_matrix_input_list.append(temp_list)
            #  gen_cartesian_input_list.append(temp_list)
            gen_cartesian_input_list.append(temp_cartesian)

    # print('gen_cartesian_input_list in gen_nrgrid_positions',gen_cartesian_input_list)
    # generate z-matrix
    return gen_cartesian_input_list


def gen_nrgrid_positions_v2(atm, link_list, value_list, rgrid_list, nth_dof):
    # e.g.,  [[], [0.49999996680528824]], 0
    # or #  ['C', 'H', 'N'] link_list [[], ['1'], ['1', '2']] value_list [[], ['1.064'], ['1.156', '180.']], 0

    print("gen_nrgrid_positions_v2 input", "rgrid_list", rgrid_list, "nth_dof", nth_dof)

    gen_cartesian_input_list = []

    # gathering z-matrix
    for i, rgrid in enumerate(rgrid_list):

        temp_value_list = copy.deepcopy(value_list)
        for j, value in enumerate(temp_value_list):
            if j == 0:
                continue
            else:
                # for (k,rsp_entry) in enumerate(output_rsp):
                # provide temp_value from rsp lists
                # N.B., so far only bond distance, no angle
                if nth_dof + 1 == j:
                    temp_value_list[j][0] = rgrid * bohr2Ang

        # temp_list.append(temp_value_list)
        # print('temp_list', temp_list)
        # print('atm',atm, 'link_list', link_list, 'temp_value_list', temp_value_list)
        print("temp_value_list", temp_value_list)
        temp_cartesian = gen_poly_cartesian(atm, link_list, temp_value_list)

        # gen_z_matrix_input_list.append(temp_list)
        #  gen_cartesian_input_list.append(temp_list)
        gen_cartesian_input_list.append(temp_cartesian)

    # print('gen_cartesian_input_list in gen_nrgrid_positions',gen_cartesian_input_list)
    # generate z-matrix
    return gen_cartesian_input_list


# gen_poly_cartesian input
# ['C', 'H', 'N']
# link_list [[], ['1'], ['1', '2']]
# value_list [[], ['1.064'], ['1.156', '180.']]
# output: z-matrix
def gen_poly_cartesian(atm, link_list, value_list):
    # print('gen_poly_cartesian input',atm)
    # print('link_list',link_list)
    # print('value_list', value_list)
    # assemble internal coordinate
    z_matrix_list = []
    for i in range(len(atm)):
        temp_line = atm[i]
        # print('lll',link_list, len(link_list[i]))
        for j, link in enumerate(link_list[i]):
            #     print('!!', j,link)
            temp_line = temp_line + " " + str(link_list[i][j])
            temp_line = temp_line + " " + str(value_list[i][j])

        #          print(link_list[j])
        #    print(temp_line)
        z_matrix_list.append(temp_line)

    # z_matrix_list.insert(0, '""" ')
    # z_matrix_list.insert(len(z_matrix_list), '""" ')
    z_matrix = "\n".join(z_matrix_list)

    # print('z_matrix',z_matrix)
    hcn = psi4.geometry(z_matrix)
    full_geometry = hcn.full_geometry()

    # key line, is this necessary?
    full_geometry_mat = np.array(full_geometry)
    # full_geometry_mat = full_geometry
    # default is in Bohr
    # print('full_geometry_mat Bohr',full_geometry_mat)
    # print(full_geometry_mat.shape)
    # print(full_geometry_mat.shape[0])

    for i in range(full_geometry_mat.shape[0]):
        for j in range(full_geometry_mat.shape[1]):
            full_geometry_mat[i][j] = (
                full_geometry_mat[i][j] * psi4.constants.bohr2angstroms
            )

    # full_geometry_mat = np.array(full_geometry_mat)
    # full_geometry_mat = list(full_geometry_mat)
    # for i in range(len(full_geometry_mat)):
    #  full_geometry_mat[i] = list(full_geometry_mat[i])

    # print('full_geometry_mat ang', full_geometry_mat, type(full_geometry_mat))
    # full_geometry_mat = adjust_coordinates(full_geometry_mat)
    return full_geometry_mat
    ##print(print_out)
    # for i in range(len(print_out)):
    #    print(i, print_out[i] )
    # hcn = psi4.geometry("""
    # C
    # H 1 1.064
    # N 1 1.156 2 180.
    # """)

    # hcn.print_out()


# sys.exit("test")


def adjust_one_atom(input_coordinates):
    # e.g., input_coordinates [[0.0, 0.0, 0.2797334848324073], [0.0, 0.0, 0.8115565815557572], [0.0, 0.0, -0.2981280292192328], [0.0, 0.0, 0.5793102421478376], [0.0, 0.0, 1.1111333388711873], [0.0, 0.0, -0.5764127859554427], [0.0, 0.0, 0.25987046779964984], [0.0, 0.0, 1.32404583845702], [0.0, 0.0, -0.31799104625199015], [0.0, 0.0, 0.5594472251150802], [0.0, 0.0, 1.6236225957724502], [0.0, 0.0, -0.5962758029882]]

    output_coordinates = copy.deepcopy(input_coordinates)
    # print(len(input_list))

    dist = np.zeros((len(input_coordinates), len(input_coordinates)))
    # print(dist.shape)

    flag_adjust = 0
    threshd_adjust = 0.1

    for i, _ in enumerate(input_coordinates):
        for j, _ in enumerate(input_coordinates):

            dist[i, j] = np.sqrt(
                (input_coordinates[i][0] - input_coordinates[j][0]) ** 2
                + (input_coordinates[i][1] - input_coordinates[j][1]) ** 2
                + (input_coordinates[i][2] - input_coordinates[j][2]) ** 2
            )
            # print(i, j, dist[i,j], type(dist[i,j]) )

            if i < j and dist[i, j] < threshd_adjust:

                flag_adjust = 1

                print(i, j, dist[i, j])
                print("before", output_coordinates[i][2], output_coordinates[j][2])
                # N.B. assume only z axis has non-zero values
                # so far only for linear molecule or non-linear but only for z axis
                if input_coordinates[i][2] < input_coordinates[j][2]:
                    # print('case a')
                    # elements_j[2] = copy.deepcopy(elements_j[2] + 0.1)
                    output_coordinates[j][2] = output_coordinates[j][2] + 0.1

                else:
                    # print('case b')
                    output_coordinates[j][2] = output_coordinates[j][2] - 0.1
                #    elements_j[2] = copy.deepcopy(elements_j[2] - 0.1)

                print("after", output_coordinates[i][2], output_coordinates[j][2])

                dist[i, j] = np.sqrt(
                    (output_coordinates[i][0] - output_coordinates[j][0]) ** 2
                    + (output_coordinates[i][1] - output_coordinates[j][1]) ** 2
                    + (output_coordinates[i][2] - output_coordinates[j][2]) ** 2
                )

                print("after", i, j, dist[i, j])

            if flag_adjust == 1:
                break

        if flag_adjust == 1:
            break

    #'output coordinate', output_coordinates)

    return output_coordinates


def print_dist_matrix(input_coordinates):

    print("print_dist_matrix")
    dist = np.zeros((len(input_coordinates), len(input_coordinates)))
    dist_off_diagonal_list = []
    print(dist.shape)

    for i, line_i in enumerate(input_coordinates):
        for j, line_j in enumerate(input_coordinates):
            dist[i, j] = np.sqrt(
                (input_coordinates[i][0] - input_coordinates[j][0]) ** 2
                + (input_coordinates[i][1] - input_coordinates[j][1]) ** 2
                + (input_coordinates[i][2] - input_coordinates[j][2]) ** 2
            )

            if i != j:
                dist_off_diagonal_list.append(dist[i, j])

    #        print('after',i,j,dist[i,j])
    print(dist)
    print("min", min(dist_off_diagonal_list))
    return min(dist_off_diagonal_list)


# def gen_distant_matrix(input_str):
def adjust_coordinates(input_list):
    print("adjust_coordinates input------------------", input_list)
    # convert into format from
    # [ [[x,y,z (atom1)], [x,y,z(atom2)]], [ molecule2], ...
    # to
    # [ [x,y,z (atom1)], [x,y,z(atom2)]
    # print(input_list, type(input_list))

    # print(type(input_list))
    # input_list = list(input_list)
    # input_list = list(input_list)
    # print(type(input_list))
    # print(input_list)
    # for i in range(len(input_list)):
    #    input_list[i] = list(input_list[i])
    # print(type(input_list),type(input_list[0]))
    # print(input_list)

    n_rsp = len(input_list)
    n_dof = len(input_list[0])

    print("n_rsp", n_rsp, "n_dof", n_dof)

    input_coordinates = []
    for i, line_i in enumerate(input_list):
        #     #line_i = list(line_i)
        #
        #     print('loop-1', i, line_i, type(line_i))
        for j, line_j in enumerate(line_i):
            #          print('loop-2', i, j, line_j, type(line_j),list(line_j), type(list(line_j)) )
            input_coordinates.append(line_j)
    # input_coordinates = copy.deepcopy(input_list)

    print("input_coordinates", input_coordinates)
    print("distance matrix!!!!!!!!!!!!!!!!!!!")
    print_dist_matrix(input_coordinates)

    # temp_coordinates = copy.deepcopy(input_coordinates)
    # print(len(input_list))
    # dist = np.zeros((len(input_coordinates),len(input_coordinates)))
    # print(dist.shape)
    # print(dist)
    # print('temp_coordinates', temp_coordinates)
    # temp_coordinates = adjust_one_atom(input_coordinates)
    iter = 16
    threshd = 0.05
    temp_coordinates = adjust_one_atom(input_coordinates)

    print("min2", print_dist_matrix(temp_coordinates))

    for i in range(iter):
        temp_coordinates = adjust_one_atom(temp_coordinates)
        if print_dist_matrix(temp_coordinates) > threshd:
            print("converged")
            break
    # temp_coordinates2 = adjust_one_atom(temp_coordinates)
    # temp_coordinates3 = adjust_one_atom(temp_coordinates2)

    print("temp_coordinates", temp_coordinates)

    # convert back
    output_list = copy.deepcopy(temp_coordinates)
    # output_list = []
    # for (i,line) in enumerate(temp_coordinates):
    #     if i%3 == 0:
    #         temp_list = []
    #         temp_list.append(line)
    #     elif i%3==1:
    #         temp_list.append(line)
    #     else:
    #         temp_list.append(line)
    #         output_list.append(temp_list)

    # need to be assembled together
    # input
    #  [array([ 0.        ,  0.        , -0.10048341]), array([0.        , 0.        , 0.69951823]), array([ 0.        ,  0.        , -0.20096683]), array([0.        , 0.        , 1.39903645]), array([ 0.        ,  0.        , -0.35169195]), array([0.        , 0.        , 2.44831379]), array([ 0.        ,  0.        , -0.50241707]), array([0.        , 0.        , 3.49759113])]
    # output  [array([[ 0.        ,  0.        , -0.10048341],
    #    [ 0.        ,  0.        ,  0.69951823]]), array([[ 0.        ,  0.        , -0.20096683],
    #    [ 0.        ,  0.        ,  1.39903645]]), array([[ 0.        ,  0.        , -0.35169195],
    #    [ 0.        ,  0.        ,  2.44831379]]), array([[ 0.        ,  0.        , -0.50241707],
    #    [ 0.        ,  0.        ,  3.49759113]])]

    print("output_list", output_list)
    output_list_2 = []
    m = 0
    for i in range(n_rsp):
        coord_temp = []
        for j in range(n_dof):
            coord_temp.append(list(output_list[m]))
            m = m + 1

        coord_temp = np.array(coord_temp)
        output_list_2.append(coord_temp)

    print("output_list_2", output_list_2)
    return output_list_2


# sys.exit()


def gen_grid_dist():
    print("test")


def gen_grid_angle():
    print("test")


def gen_grid_dihedral():
    print("test")
