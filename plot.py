import csv
import os
import math

start_value = 15

def clear_file(file):
    """
    removes all `?` from the beginning of an file
    :param file:
    :return:
    """
    with open(file) as f:
        data = []
        while True:
            s = f.readline()
            if not s:
                break

            if len(s) == 0:
                continue

            data.append(s.replace("?", ""))
    return data


def number_of_relations_found(file):
    data = clear_file(file)
    reader = csv.reader(data, delimiter=',', quotechar='|')

    l_data = []
    d_data = []
    for row in reader:
        m         = float(row[4])
        relations = float(row[8])/float(row[10]) # align to benches
        if row[0] == "l":
            l_data.append(relations/m)
        else:
            d_data.append(relations/m)

    print(l_data)
    print(d_data)

    return sum(l_data)/len(l_data), sum(d_data)/len(d_data)


def read_csv(file, aslog=True, index=11, diff=-1):
    """
    reads the csv file from mem
    :param data: =
        # first is legendre runtime,
        ({
            "m":   [0.1, 0.2, 0.3, ..., 0.9],
            "m+1": [0.1, 0.2, 0.3, ..., 0.9],
            ...
        },
        # this is the dlog runtime
        {
            "m":   [0.1, 0.2, 0.3, ..., 0.9],
            "m+1": [0.1, 0.2, 0.3, ..., 0.9],
            ...
        })
    :param aslog: if the data
    :param index: index of the time index
    :param diff: divide the time number
    :return:
    """
    data = clear_file(file)

    i_m = 4
    i_s = 3
    i_t = 2
    i_time = index
    diff_ = 1.
    l_ret_data = {}
    d_ret_data = {}

    reader = csv.reader(data, delimiter=',', quotechar='|')
    for row in reader:
        m = row[i_m]
        if m not in l_ret_data:
            l_ret_data[m] = []
            d_ret_data[m] = []

        if diff != -1:
            diff_ = float(row[diff])

        if aslog:
            if row[0] == "l":
                l_ret_data[m].append(math.log(float(row[i_time])/diff_, 2))
            else:
                d_ret_data[m].append(math.log(float(row[i_time])/diff_, 2))
        else:
            if row[0] == "l":
                l_ret_data[m].append(float(row[i_time])/diff_)
            else:
                d_ret_data[m].append(float(row[i_time])/diff_)

    return l_ret_data, d_ret_data


def read_csv_no_m(file, aslog=True, index=11, diff=-1, mult=-1):
    """
    reads the csv file from mem
    :param data: =
        # first is legendre runtime,
        ({
            "m":   [0.1, 0.2, 0.3, ..., 0.9],
            "m+1": [0.1, 0.2, 0.3, ..., 0.9],
            ...
        },
        # this is the dlog runtime
        {
            "m":   [0.1, 0.2, 0.3, ..., 0.9],
            "m+1": [0.1, 0.2, 0.3, ..., 0.9],
            ...
        })
    :param aslog: if the data
    :param index: index of the time index
    :param diff: divide the time number
    :return:
    """
    data = clear_file(file)

    i_m = 4
    i_s = 3
    i_t = 2
    i_time = index
    diff_ = 1.
    mult_ = 1.
    l_ret_data = {}
    d_ret_data = {}

    l_ret_data[0] = []
    d_ret_data[0] = []

    reader = csv.reader(data, delimiter=',', quotechar='|')
    for row in reader:
        if diff != -1:
            diff_ = float(row[diff])

        if mult != -1:
            mult_ = float(row[mult])

        if aslog:
            if row[0] == "l":
                l_ret_data[0].append(math.log(float(row[i_time])*mult_/diff_, 2))
            else:
                d_ret_data[0].append(math.log(float(row[i_time])*mult_/diff_, 2))
        else:
            if row[0] == "l":
                l_ret_data[0].append(float(row[i_time])/diff_)
            else:
                d_ret_data[0].append(float(row[i_time])/diff_)

    return l_ret_data, d_ret_data


def linear_interpolate(data):
    from scipy.stats import linregress
    from copy import deepcopy

    retdata = deepcopy(data)

    function_data = []
    for tup in range(2):
        d = data[tup]
        for m, v in d.items():
            l = len(v)
            X = list(range(l))
            Y = [float(vv) for vv in v]
            b, a, r, p, std = linregress(X, Y)
            # flingress = a + b*x
            function_data.append([a, b])
            # print(b, a, r, p, std, "f =", str(a) + " + " + str(b) + "*x")
            retdata[tup][m].clear()
            retdata[tup][m] = [(a + (b * x)) for x in X]

    return retdata


def write_tex(out_name, data, single_m=False):
    """
    writes `data` to `out_name` in a text readable format.
        If `single_m` is set this function will not generate for each different
        `m` a new output file.
    :param out_name
    :param data:
    :param single_m
    :return:
    """
    split = os.path.split(out_name)

    for tup, fname in enumerate([split[0] + "/l_m" + "_" + split[1],
                                 split[0] + "/d_m" + "_" + split[1]]):
        for m, _ in data[tup].items():
            with open(fname, "w") as f:
                f.write("W T\n")
                i = start_value
                for v in data[tup][m]:
                    f.write(str(i) + " " + str(v) + "\n")
                    i += 1

# data_prep_m1_s13_t13 = read_csv("data/prep_m1_s13_t13.log")
# data_mult_prep_m112_s13 = read_csv_no_m("data/prep_m112_s13_t13.log", True, 11, -1, 4)
# data_mult_prep_m112_s13 = read_csv_no_m("data/prep_m112_s13_t13.log")
# data_mult_prep_m112_s13 = read_csv_no_m("data/prep_m112_s13_t13_iters100.log")
# data_mult_prep_m16_s13 = read_csv_no_m("data/multprep_m16_s13.log")
# data_mult_prep_m16_s13 = read_csv_no_m("data/multprep_m16_s13_v2.log",  True, 11, -1, 3)
# data_mult_prep_m16_s13 = read_csv_no_m("data/multprep_m16_s13_iters100.log")
# data_mult_prep_m16_s13 = read_csv_no_m("data/multprep_m16_s13_iters100.log", True, 11, -1, 4)
# data_noprep_m112 = read_csv_no_m("data/noprep_m112_s13_t13.log", True, 7, 6)
# data_noprep_m13 = read_csv_no_m("data/noprep_13.log", True, 7, 6)
# data_noprep_m13 = read_csv_no_m("data/noprep_13_v2.log", True, 7, 6)
# data_noprep_m13 = read_csv_no_m("data/noprep_13_iters1_fixxed_dlog2.log", True, 7, 6)

# #data_prep_m1_s13_t13 = read_csv("data/prep_m1_s13_t13_iters100.log")
# #data_prep_m112_s13_t13 = read_csv_no_m("data/prep_m112_s13_t13_iters100.log")
# #data_noprep_m112_s13_t13 = read_csv_no_m("data/noprep_m112_s13_t13_iters100.log", True, 7, 6)
#
# write_tex("data/.out/prep_m1_s13", data_prep_m1_s13_t13)
# write_tex("data/.out/mult_prep_m112_s13", data_mult_prep_m112_s13, True)
# write_tex("data/.out/mult_prep_m16_s13", data_mult_prep_m16_s13, True)
# write_tex("data/.out/noprep_m112_s13", data_noprep_m112_s13, True)
# write_tex("data/.out/noprep_m13", data_noprep_m13, True)
#
# write_tex("data/.out/interpolate_prep_m1_s13", linear_interpolate(data_prep_m1_s13_t13))
# write_tex("data/.out/interpolate_mult_prep_m112_s13", linear_interpolate(data_mult_prep_m112_s13), True)
# write_tex("data/.out/interpolate_mult_prep_m16_s13", linear_interpolate(data_mult_prep_m16_s13), True)
# #write_tex("data/.out/interpolate_noprep_m112_s13", linear_interpolate(data_noprep_m112), True)
# write_tex("data/.out/interpolate_noprep_m13", linear_interpolate(data_noprep_m13), True)

# print("%relations found: ", number_of_relations_found("data/noprep_13_v2.log"))
