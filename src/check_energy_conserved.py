import sys

Etot = []
FILENAME_0 = sys.argv[1]
threshold = 1.0e-3
with open(FILENAME_0) as f_0:
    LINES_0 = f_0.readlines()
    for n, line in enumerate(LINES_0):
        if n in [0]:
            continue

        if n == 1:
            terms = line.split()
            # print(terms)
            for m, term in enumerate(terms):
                if term == "Htot":
                    loc = m

        # print(loc)
        else:
            terms = line.split()
            Etot.append(float(terms[loc]))

    # print(Etot)
    Emax = max(Etot)
    Emin = min(Etot)
    # print(Emax,Emin)
    if (Emax - Emin) > threshold:
        print("energy is not conserved")
    else:
        print("energy is conserved")
