import csv
import matplotlib.pyplot as plt
import matplotlib
from itertools import cycle

matplotlib.rcParams['figure.figsize'] = 26, 10
cycol = cycle('bgrcmk')

Results = {}
with open("result.log") as log:
    current = ""
    for line in log:
        line = line.strip()
        if not line:
            continue
        parts = line.split(" ")
        if parts[2] == "cannot" or parts[2].isdigit():
            if parts[2].isdigit():
                Results[current].append({
                    "matrix": parts[0],
                    "initial nnz": int(parts[1]),
                    "result nnz": int(parts[2]),
                    "total ms": int(parts[3]),
                })
            else:
                Results[current].append({
                    "matrix": parts[0],
                    "initial nnz": int(parts[1]),
                    "result nnz": -1,
                    "total ms": -1,
                })
        else:
            current = line
            Results[current] = []

cur_itr = 0
fig, ax = plt.subplots(4, 2)
for key, dict_arr in Results.items():
    if key[:7] == "Reorder":
        vals_for_ini = {}
        for dot in dict_arr:
            if dot["total ms"] != -1:
                vals_for_ini[dot["total ms"]] = dot["result nnz"]

        vals_for_ini = dict(sorted(vals_for_ini.items()))
        ax[cur_itr, 0].plot(vals_for_ini.keys(), vals_for_ini.values(),
                         label="Result nnz "+key, color=next(cycol))
        ax[cur_itr, 0].set_xticks(range(0, 3400001, 200000))
        ax[cur_itr, 0].set_yticks(range(0, 400000001, 50000000))
        cur_itr += 1

cur_itr = 0
for key, dict_arr in Results.items():
    if key[:7] == "Reorder":
        vals_for_ini = {}
        for dot in dict_arr:
            if dot["total ms"] != -1:
                vals_for_ini[dot["total ms"]] = dot["initial nnz"]

        vals_for_ini = dict(sorted(vals_for_ini.items()))
        ax[cur_itr, 1].plot(vals_for_ini.keys(), vals_for_ini.values(),
                         label="Initial nnz "+key, color=next(cycol))
        ax[cur_itr, 1].set_xticks(range(0, 3400001, 200000))
        ax[cur_itr, 1].set_yticks(range(0, 16000001, 2000000))
        cur_itr += 1

fig.tight_layout()
fig.legend(loc='lower center', ncol=8)
#plt.show()
plt.savefig('result.png')


with open("result.csv", "w", newline="") as csvfile:
    fieldnames = ["matrix", "initial nnz", "result nnz", "total ms"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    for key, dict_arr in Results.items():
        csvfile.write(key+",\n")
        writer.writeheader()
        writer.writerows(dict_arr)
