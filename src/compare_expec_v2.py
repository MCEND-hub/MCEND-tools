"""
This code checks the consistency of data in two expect.t generated by MCEND.

Criteria:
    1. If two expect.t have the same number of lines
    2. If each line has the same number of entries
    3. If the absolute values of an entry in two expec.t are smaller than a given threshold,
        consistent
    4. If at least one entry of the absolute value greater than the given thresold,
       compare the difference of the two entries, if the difference below the given threshold,
       consistent

Usage: (assume two expec files below exist)
python compare_expec_v2.py expec.t-lih-cs ./calc_testing/expec.t-lih-cs-fc-ref

Cong Wang
Oct-10-2023
"""

import sys


FILENAME_0 = sys.argv[1]
FILENAME_1 = sys.argv[2]

THRESHOLD_1 = 1.0e-5

FLAG_G = 0

with open(FILENAME_0, encoding="utf-8") as f_0:
    LINES_0 = f_0.readlines()

with open(FILENAME_1, encoding="utf-8") as f_1:
    LINES_1 = f_1.readlines()

    if len(LINES_1) != len(LINES_0):
        # rule 1
        print("inconsistent number of lines in two files")
        FLAG_G = 1

    for n, line in enumerate(LINES_1):
        # skip title lines
        if n in [0, 1]:
            continue

        data_0_ = LINES_0[n].split(" ")
        data_0_ = list(filter(None, data_0_))
        data_0 = [float(i) for i in data_0_]

        data_1_ = LINES_1[n].split(" ")
        data_1_ = list(filter(None, data_1_))
        data_1 = [float(i) for i in data_1_]

        if len(data_0) != len(data_1):
            # rule 2
            print("inconsistent number of entries in line", n)
            FLAG_G = 1

        for m, datum in enumerate(data_0):
            flag = 0

            if abs(data_0[m]) < THRESHOLD_1 and abs(data_1[m]) < THRESHOLD_1:
                # print(' skip small entries', n, LINES_0[n])
                continue

            diff = abs(data_0[m]) - abs(data_1[m])
            diff = abs(diff)

            if diff >= THRESHOLD_1:
                FLAG_G = 1
                print("WARNING", n, data_0[m], data_1[m])

if FLAG_G == 0:
    print("Consistent")
