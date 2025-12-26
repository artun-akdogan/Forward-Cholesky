import csv
import re

input_file_mat = "result_mat.log"
input_file_oct = "result_oct.log"
input_file_my = "result_my.log"
output_file = "output.csv"


"""{
    "filename": None,
    "type": None,
    "elapsed_time": None,
    "begin_nnz": None,
    "result_nnz": None,
    "mine_nnz": None,
    "frobenius": None,
    "exit_status": None
}"""


def parse_mattype(input_file):

    data = {}

    # Read input file
    with open(input_file, "r", encoding="utf-8") as file:
        content = file.read().split("-----------\n")

    lookup = set()

    for item in content:
        if item.strip() == "":
            continue
        rows = item.strip().split("\n")
        res = {
            "filename": None,
            "type": None,
            "elapsed_time": None,
            "begin_nnz": None,
            "result_nnz": None,
            "mine_nnz": None,
            "frobenius": None,
            "exit_status": None
        }
        for itr in range(len(rows)):
            if rows[itr].strip().startswith("matrices/"):
                res["filename"] = rows[itr].strip().partition("/")[2]
            elif rows[itr].strip().startswith("Reorder"):
                if res["filename"] in lookup:
                    res["type"] = "Reorder Sequential"
                else:
                    lookup.add(res["filename"])
                    res["type"] = "Reorder Parallel"
            elif rows[itr].strip().startswith("Original"):
                if res["filename"] in lookup:
                    res["type"] = "Original Sequential"
                else:
                    res["type"] = "Original Parallel"
            elif rows[itr].strip().startswith("Elapsed time is "):
                res["elapsed_time"] = float(rows[itr].strip().split()[3])
            elif rows[itr].strip().startswith("Begin nnz"):
                try:
                    res["begin_nnz"] = int(rows[itr+1].strip())
                except ValueError:
                    res["begin_nnz"] = float(rows[itr+1].strip())
            elif rows[itr].strip().startswith("Result nnz"):
                try:
                    res["result_nnz"] = int(rows[itr+1].strip())
                except ValueError:
                    res["result_nnz"] = float(rows[itr+1].strip())
            elif rows[itr].strip().startswith("Mine nnz"):
                try:
                    res["mine_nnz"] = int(rows[itr+1].strip())
                except ValueError:
                    res["mine_nnz"] = float(rows[itr+1].strip())
            elif rows[itr].strip().startswith("Frobenius Scalar Product"):
                res["frobenius"] = float(rows[itr+1].strip())
            elif rows[itr].strip().startswith("EXIT"):
                res["exit_status"] = int(
                    rows[itr].strip().rpartition(" ")[2][:-1])
        if res["filename"] not in data:
            data[res["filename"]] = {"Reorder Sequential": {}, "Reorder Parallel": {
            }, "Original Sequential": {}, "Original Parallel": {}}
        data[res["filename"]][res["type"]] = res

    return data


def parse_my(input_file):

    data = {}

    # Read input file
    with open(input_file, "r", encoding="utf-8") as file:
        content = file.read().split("-----------\n")

    for item in content:
        if item.strip() == "":
            continue
        rows = item.strip().split("\n")
        res = {
            "filename": None,
            "type": None,
            "elapsed_time": None,
            "mine_nnz": None,
            "exit_status": None
        }
        for itr in range(len(rows)):
            if rows[itr].strip().startswith("Start for matrices/"):
                res["filename"] = rows[itr].strip().partition("/")[2]
            elif rows[itr].strip().startswith("Reorder Sequential"):
                res["type"] = "Reorder Sequential"
            elif rows[itr].strip().startswith("Reorder Parallel"):
                res["type"] = "Reorder Parallel"
            elif rows[itr].strip().startswith("Original Sequential"):
                res["type"] = "Original Sequential"
            elif rows[itr].strip().startswith("Original Parallel"):
                res["type"] = "Original Parallel"
            elif rows[itr].strip().startswith("Total time except disk io:"):
                res["elapsed_time"] = float(
                    rows[itr].strip().split()[5]) / 1000.0
            elif rows[itr].strip().startswith("Total nonzero: "):
                res["mine_nnz"] = int(rows[itr].strip().split()[2])
            elif rows[itr].strip().startswith("EXIT"):
                res["exit_status"] = int(
                    rows[itr].strip().rpartition(" ")[2][:-1])

        if res["filename"] not in data:
            data[res["filename"]] = {"Reorder Sequential": {}, "Reorder Parallel": {
            }, "Original Sequential": {}, "Original Parallel": {}}
        data[res["filename"]][res["type"]] = res

    return data


data_mat = parse_mattype(input_file_mat)
data_oct = parse_mattype(input_file_oct)
data_my = parse_my(input_file_my)

csv_data = {"Original": [], "Reorder": []}
csv_ctrl = {"Original": [], "Reorder": []}


def from_stat(exit_code):
    if exit_code == 124:
        return "*t"
    if exit_code == 137:
        return "*m"
    return "*u"


for name in data_my:
    for typ in ["Original", "Reorder"]:
        result = {
            "Matrix Name": name,
            "Column Size": None,
            "Begin NNZ": None,
            "Sequential Time": None,
            "Parallel Time": None,
            "Matlab Time": None,
            "Octave Time": None,
            #"Algo. Success": None
        }
        res_ctrl = {
            "Matrix Name": name,
            "Begin NNZ": None,
            "Algo. NNZ": None,
            "Matlab NNZ": None,
            #"NNZ Ratio": None,
            "Forward/Matlab Norm": None,
            "Matlab/Octave Norm": None,
        }
        nnz = 0
        print(name, typ)
        with open(f"matrices/{name}", "r") as fil:
            lines = fil.readlines()
            read_once = False
            for line in lines:
                if line.strip().startswith("%"):
                    continue
                if not read_once:
                    result["Column Size"] = int(line.split()[0])
                    read_once = True
                    continue
                try:
                    if float(line.split(" ")[2])!= 0.0:
                        nnz+=1
                except ValueError:
                    pass
        
        res_ctrl["Begin NNZ"] = nnz
        result["Begin NNZ"] = nnz
        if data_my[name][f"{typ} Sequential"]["exit_status"] == 0:
            result["Sequential Time"] = data_my[name][f"{typ} Sequential"]["elapsed_time"]
            res_ctrl["Algo. NNZ"] = data_my[name][f"{typ} Sequential"]["mine_nnz"]
        else:
            result["Sequential Time"] = from_stat(
                data_my[name][f"{typ} Sequential"]["exit_status"])
        if data_my[name][f"{typ} Parallel"]["exit_status"] == 0:
            result["Parallel Time"] = data_my[name][f"{typ} Parallel"]["elapsed_time"]
            if res_ctrl["Algo. NNZ"] is None:
                res_ctrl["Algo. NNZ"] = data_my[name][f"{typ} Parallel"]["mine_nnz"]
            else:
                assert res_ctrl["Algo. NNZ"] == data_my[name][f"{typ} Parallel"]["mine_nnz"]
        else:
            result["Parallel Time"] = from_stat(
                data_my[name][f"{typ} Parallel"]["exit_status"])
            res_ctrl["Algo. NNZ"] = from_stat(
                data_my[name][f"{typ} Parallel"]["exit_status"])

        if data_mat[name][f"{typ} Parallel"]["exit_status"] == 0:
            result["Matlab Time"] = data_mat[name][f"{typ} Parallel"]["elapsed_time"]
            if res_ctrl["Begin NNZ"] is None:
                res_ctrl["Begin NNZ"] = data_mat[name][f"{typ} Parallel"]["begin_nnz"]
            else:
                assert res_ctrl["Begin NNZ"] == (data_mat[name][f"{typ} Parallel"]["begin_nnz"]+result["Column Size"])//2

            if res_ctrl["Algo. NNZ"] is None:
                pass
                #res_ctrl["Algo. NNZ"] = data_mat[name][f"{typ} Parallel"]["mine_nnz"]
            else:
                pass
            """
                try:
                    assert res_ctrl["Algo. NNZ"] == data_mat[name][f"{typ} Parallel"]["mine_nnz"]
                except Exception as exc:
                    print(res_ctrl["Algo. NNZ"], data_mat[name][f"{typ} Parallel"]["mine_nnz"])
                    raise exc
            """
            res_ctrl["Matlab NNZ"] = data_mat[name][f"{typ} Parallel"]["result_nnz"]
            if res_ctrl["Algo. NNZ"] is not None:
                pass
                #result["NNZ Ratio"] = result["Algo. NNZ"] / result["Matlab NNZ"]
            if data_mat[name][f"{typ} Parallel"]["frobenius"] is not None:
                res_ctrl["Forward/Matlab Norm"] = "{:.3e}".format(data_mat[name][f"{typ} Parallel"]["frobenius"])
        else:
            result["Matlab Time"] = from_stat(
                data_mat[name][f"{typ} Parallel"]["exit_status"])
            res_ctrl["Matlab NNZ"] = from_stat(
                data_mat[name][f"{typ} Parallel"]["exit_status"])
        if data_oct[name][f"{typ} Parallel"]["exit_status"] == 0:
            result["Octave Time"] = data_oct[name][f"{typ} Parallel"]["elapsed_time"]

            if data_oct[name][f"{typ} Parallel"]["frobenius"] is not None:
                res_ctrl["Matlab/Octave Norm"] = "{:.3e}".format(data_oct[name][f"{typ} Parallel"]["frobenius"])
        else:
            result["Octave Time"] = from_stat(
                data_oct[name][f"{typ} Parallel"]["exit_status"])

        if (isinstance(result["Parallel Time"], float) and not isinstance(result["Octave Time"], float)) or (isinstance(result["Parallel Time"], float) and result["Parallel Time"] < result["Octave Time"]):
            pass
            #result["Algo. Success"] = True
        csv_data[typ].append(result)
        csv_ctrl[typ].append(res_ctrl)


# Write to CSV
with open("output_original.csv", "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_data["Original"][0].keys())
    writer.writeheader()
    writer.writerows(csv_data["Original"])
print(f"Data extracted and saved to output_original.csv")

with open("output_reorder.csv", "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_data["Reorder"][0].keys())
    writer.writeheader()
    writer.writerows(csv_data["Reorder"])
print(f"Data extracted and saved to output_reorder.csv")

with open("ctrl_original.csv", "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_ctrl["Original"][0].keys())
    writer.writeheader()
    writer.writerows(csv_ctrl["Original"])
print(f"Data extracted and saved to ctrl_original.csv")

with open("ctrl_reorder.csv", "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_ctrl["Reorder"][0].keys())
    writer.writeheader()
    writer.writerows(csv_ctrl["Reorder"])

print(f"Data extracted and saved to ctrl_reorder.csv")
