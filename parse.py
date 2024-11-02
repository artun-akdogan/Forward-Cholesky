import csv
import matplotlib.pyplot as plt

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

vals_for_ini = {}
vals_for_res = {}
for key, dict_arr in Results.items():
    for dot in dict_arr:
        if dot["total ms"] != -1:
            vals_for_ini[dot["total ms"]] = dot["initial nnz"]
            vals_for_res[dot["total ms"]] = dot["result nnz"]

fig, ax = plt.subplots(1, 2)
vals_for_ini = dict(sorted(vals_for_ini.items()))
ax[0].plot(vals_for_ini.keys(), vals_for_ini.values(), label="initial nnz", color="green")
vals_for_res = dict(sorted(vals_for_res.items()))
ax[1].plot(vals_for_res.keys(), vals_for_res.values(), label="result nnz", color="blue")
fig.legend()
plt.savefig('result.png')


with open("result.csv", "w", newline="") as csvfile:
    fieldnames = ["matrix", "initial nnz", "result nnz", "total ms"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    for key, dict_arr in Results.items():
        csvfile.write(key+",\n")
        writer.writeheader()
        writer.writerows(dict_arr)
