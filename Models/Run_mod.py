from argparse import ArgumentParser
import subprocess
import pandas as pd
import re
import time

def parse_slim(slim_string, drive_number):
    """
    Parse SLiM output and return all recorded information as lists.
    """
    gen, total, female, rate_homo_fe, rate_hete_fe = [], [], [], [], []
    wt = [[] for _ in range(drive_number)]
    dr = [[] for _ in range(drive_number)]
    r1 = [[] for _ in range(drive_number)]
    r2 = [[] for _ in range(drive_number)]
    rate, new_infections, new_1, new_2, new_3, new_4, pop_immune = [], [], [], [], [], [], []
    native_new_infections, intro_human = [], []
    ann, vnn = [], []
    suppressed, gen_suppressed = 0, 10000
    mark1 = mark2 = mark3 = None
    current_ann, current_vnn = 0, 0

    for line in slim_string.split('\n'):
        spaced = line.split()

        # Suppression info
        if line.startswith("SUPPRESSED_WITHOUT_CHASING:: "):
            suppressed, gen_suppressed = 1, int(spaced[1])
        elif line.startswith("SUPPRESSED_WITH_CHASING:: "):
            suppressed, gen_suppressed = 2, int(spaced[1])
        elif line.startswith("SUPPRESSED:: "):
            suppressed, gen_suppressed = 1, int(spaced[1])

        elif line.startswith("SPATIAL_GEN:: "):
            current_ann = float(spaced[3]) if spaced[3] != "NULL" else 0
            current_vnn = float(spaced[5]) if spaced[5] != "NULL" else 0

        # MOS output
        elif line.startswith("MOS:: "):
            gen.append(int(spaced[1]))
            total.append(int(spaced[2]))
            female.append(int(spaced[3]))
            rate_homo_fe.append(float(spaced[4]))
            rate_hete_fe.append(float(spaced[5]))
            mark1, mark2, mark3 = spaced[6], spaced[7], spaced[8]
            ann.append(current_ann)
            vnn.append(current_vnn)

        # DRIVE output
        elif len(spaced) > 0 and re.match(r"^DRIVE(\d+)", spaced[0]):
            idx = int(spaced[0].split("::")[0].replace("DRIVE", ""))
            if 1 <= idx <= drive_number:
                wt[idx - 1].append(float(spaced[1]))
                dr[idx - 1].append(float(spaced[2]))
                r1[idx - 1].append(float(spaced[3]))
                r2[idx - 1].append(float(spaced[4]))

        # DENGUE output
        elif line.startswith("DENGUE:: "):
            native_new_infections.append(int(spaced[1]))
            intro_human.append(int(spaced[2]))
            rate.append(float(spaced[3]))
            new_infections.append(int(spaced[4]))
            new_1.append(int(spaced[5]))
            new_2.append(int(spaced[6]))
            new_3.append(int(spaced[7]))
            new_4.append(int(spaced[8]))
            pop_immune.append(float(spaced[10]))

    return {
        "gen": gen,
        "total": total,
        "female": female,
        "rate_homo_fe": rate_homo_fe,
        "rate_hete_fe": rate_hete_fe,
        "wt": wt,
        "dr": dr,
        "r1": r1,
        "r2": r2,
        "rate": rate,
        "new_infections": new_infections,
        "new_1": new_1,
        "new_2": new_2,
        "new_3": new_3,
        "new_4": new_4,
        "pop_immune": pop_immune,
        "intro_human": intro_human,
        "native_new_infections": native_new_infections,
        "ann": ann,
        "vnn": vnn,
        "suppressed": suppressed,
        "gen_suppressed": gen_suppressed,
        "mark1": mark1,
        "mark2": mark2,
        "mark3": mark3
    }

def summarize_slim(parsed):
    """
    Compute summary statistics from parsed SLiM output.
    """
    start_index = 15  # from week 16 onwards
    adult = parsed["total"]
    adult_fe = parsed["female"]
    rate_homo_fe = parsed["rate_homo_fe"]
    rate_hete_fe = parsed["rate_hete_fe"]
    infection = parsed["new_infections"]
    native_new_infections = parsed["native_new_infections"]
    intro_human = parsed["intro_human"]
    pop_immune = parsed["pop_immune"]

    avg_adult = sum(adult[start_index:]) / len(adult[start_index:]) if len(adult) > start_index else 0
    avg_adult_fe = sum(adult_fe[start_index:]) / len(adult_fe[start_index:]) if len(adult_fe) > start_index else 0
    avg_homo_fe = sum(rate_homo_fe[start_index:]) / len(rate_homo_fe[start_index:]) if len(rate_homo_fe) > start_index else 0
    total_infection = sum(infection[start_index:]) if len(infection) > start_index else 0
    total_native_new_infections = sum(native_new_infections[start_index:]) if len(native_new_infections) > start_index else 0
    total_intro_human = sum(intro_human[start_index:]) if len(intro_human) > start_index else 0
    immune_diff = pop_immune[-1] - pop_immune[start_index] if len(pop_immune) > start_index else 0

    return (
        parsed["suppressed"],
        parsed["gen_suppressed"],
        avg_adult, avg_adult_fe, rate_homo_fe,
        total_infection, total_intro_human, total_native_new_infections,
        parsed["mark1"], parsed["mark2"], parsed["mark3"], immune_diff
    )

def writeOutput(file, parsed, drive_number):
    headers = [
        "Generation", "Total", "Female", "Rate_homo_fe", "Rate_hete_fe",
        "Rate", "New Infections", "New 1", "New 2", "New 3", "New 4",
        "Pop Immune", "Intro Human", "Native Infections"
    ]

    drive_headers = [f"{x}_{i}" for i in range(1, drive_number + 1) for x in ("WT", "DR", "R1", "R2")]
    headers += drive_headers + ["Mark1", "Mark2", "Mark3", "ANN", "VNN"]

    df = pd.DataFrame({
        "Generation": parsed["gen"],
        "Total": parsed["total"],
        "Female": parsed["female"],
        "Rate_homo_fe": parsed["rate_homo_fe"],
        "Rate_hete_fe": parsed["rate_hete_fe"],
        "Rate": parsed["rate"],
        "New Infections": parsed["new_infections"],
        "New 1": parsed["new_1"],
        "New 2": parsed["new_2"],
        "New 3": parsed["new_3"],
        "New 4": parsed["new_4"],
        "Pop Immune": parsed["pop_immune"],
        "Intro Human": parsed["intro_human"],
        "Native New Infections": parsed["native_new_infections"]
    })
    # Add drives
    for i in range(drive_number):
        df[f"WT_{i+1}"] = parsed["wt"][i]
        df[f"DR_{i+1}"] = parsed["dr"][i]
        df[f"R1_{i+1}"] = parsed["r1"][i]
        df[f"R2_{i+1}"] = parsed["r2"][i]
    df["Mark1"] = parsed["mark1"]
    df["Mark2"] = parsed["mark2"]
    df["Mark3"] = parsed["mark3"]
    df["ANN"] = parsed["ann"]
    df["VNN"] = parsed["vnn"]

    df.to_excel(file, index=False, header=True)

def writeSummary(file, summary_tuple):
    with open(file, 'a') as f:
        f.write(','.join(map(str, summary_tuple)) + '\n')

def run_slim(cmd):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    out, err = proc.communicate()
    print(out)
    if err:
        print("Errors:", err)
    return out


def configure_slim_command_line(args_dict):
    clargs = []  # 改用列表构建，避免手动 split() 带来的歧义
    clargs.append("slim")
    source = args_dict.pop("source")

    for arg, val in args_dict.items():
        clargs.append("-d")
        if isinstance(val, bool):
            clargs.append(f"{arg}={'T' if val else 'F'}")
        elif isinstance(val, list):
            # 关键修改：移除 join 里的空格，确保 c(1,2,3,4) 是一个整体
            list_str = "c(" + ",".join(map(str, val)) + ")"
            clargs.append(f"{arg}={list_str}")
        else:
            clargs.append(f"{arg}={val}")

    clargs.append(source)
    return clargs  # 直接返回列表，不再需要 split()

def main():
    parser = ArgumentParser()
    parser.add_argument('-src', '--source', default="Mod_equ.slim", type=str,
                        help=r"SLiM file to be run.")

    parser.add_argument('-adult_females', '--NUM_ADULT_FEMALES', default=40, type=int)

    parser.add_argument('-drop', '--DROP', default=True)
    parser.add_argument('-drive_number', '--DRIVE_NUMBER', default=1, type=int)
    parser.add_argument('-drop_time', '--DROP_TIME', default=16, type=int)
    parser.add_argument('-drop_int', '--DROP_INTERVAL', default=52, type=int)

    parser.add_argument('-homo', '--HOMO_IMMUNE', default=0.95, type=float)
    parser.add_argument('-hete', '--HETE_IMMUNE', default=0.95, type=float)

    parser.add_argument('-hl', '--HAPLOLETHAL', default=False)
    parser.add_argument('-hs', '--HAPLOSUFFICIENT', default=False)

    parser.add_argument('-type', '--IMMUNE_TYPE', default=[1, 2, 3, 4])

    parser.add_argument('-drop_1', '--DROP_1', default=0.005, type=float)
    parser.add_argument('-conversion_1', '--DRIVE_CONVERSION_RATE_1', default=0.8, type=float)
    parser.add_argument('-germline_1', '--LATE_GERMLINE_RESISTANCE_CUT_RATE_1', default=0.5, type=float)
    parser.add_argument('-embryo_1', '--EMBRYO_RESISTANCE_CUT_RATE_1', default=0.1, type=float)
    parser.add_argument('-r1_1', '--R1_OCCURRENCE_RATE_1', default=0.00, type=float)
    parser.add_argument('-m_1', '--SOMATIC_FITNESS_MULTIPLIER_M_1', default=1.0, type=float)
    parser.add_argument('-f_1', '--SOMATIC_FITNESS_MULTIPLIER_F_1', default=1.0, type=float)
    parser.add_argument('-dd', '--DD_FITNESS_VALUE_1', default=0.85, type=float)
    parser.add_argument('-outfile', '--output_file', default="R1_test_mod.xlsx", type=str)
    parser.add_argument('-file', '--file1', default="Process_mod.xlsx", type=str)

    parser.add_argument('-r1_intro', '--R1_INTRODUCTION_YEAR', default=1, type=int)
    args = vars(parser.parse_args())

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    outfile = f"{args.pop('output_file')}_{timestamp}.xlsx"
    file = args.pop("file1")
    clargs = configure_slim_command_line(args)
    slim_result = run_slim(clargs)

    if ("NOT:: " not in slim_result) and ("SUPPRESSED" not in slim_result):
        print("Keywords NOT:: or SUPPRESSED not detected. Writing RE.xlsx only...")

        parsed_temp = parse_slim(slim_result, args["DRIVE_NUMBER"])
        new_row = {
            "mark1": parsed_temp["mark1"],
            "mark2": parsed_temp["mark2"],
            "mark3": parsed_temp["mark3"]
        }

        re_file = "RE.xlsx"

        try:
            # 如果文件存在，读入后追加
            df_old = pd.read_excel(re_file)
            df_new = pd.concat([df_old, pd.DataFrame([new_row])], ignore_index=True)
        except FileNotFoundError:
            # 若文件不存在，创建新文件
            df_new = pd.DataFrame([new_row])

        df_new.to_excel(re_file, index=False)

        print("Row appended to RE.xlsx. No result or summary file produced.")
        return

    parsed = parse_slim(slim_result, args["DRIVE_NUMBER"])
    summary = summarize_slim(parsed)
    writeOutput(outfile, parsed, args["DRIVE_NUMBER"])
    writeSummary(file, summary)

if __name__ == "__main__":
    main()
