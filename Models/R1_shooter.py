from argparse import ArgumentParser
import subprocess
import pandas as pd
import re
import time

def parse_slim(slim_string, drive_number):
    """
    Parse SLiM output and return all recorded information as lists.
    """
    gen, total, female, fertile_female, infertile_female = [], [], [], [], []
    wt = [[] for _ in range(drive_number)]
    dr = [[] for _ in range(drive_number)]
    r1 = [[] for _ in range(drive_number)]
    r2 = [[] for _ in range(drive_number)]
    r1_drive = [[] for _ in range(drive_number)]
    con, female_fit, r1_rate = [], [], []
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
            fertile_female.append(int(spaced[4]))
            infertile_female.append(int(spaced[5]))
            con.append(float(spaced[6]))
            female_fit.append(float(spaced[7]))
            r1_rate.append(float(spaced[8]))
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
                r1_drive[idx - 1].append(float(spaced[5]))

        # DENGUE output
        elif line.startswith("DENGUE:: "):
            native_new_infections.append(int(spaced[1]))
            intro_human.append(int(spaced[3]))
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
        "fertile_female": fertile_female,
        "infertile_female": infertile_female,
        "wt": wt,
        "dr": dr,
        "r1": r1,
        "r2": r2,
        "r1_drive": r1_drive,
        "con": con,
        "female_fit": female_fit,
        "r1_rate": r1_rate,
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
    fertile_fe = parsed["fertile_female"]
    infection = parsed["new_infections"]
    native_new_infections = parsed["native_new_infections"]
    intro_human = parsed["intro_human"]
    pop_immune = parsed["pop_immune"]

    avg_adult = sum(adult[start_index:]) / len(adult[start_index:]) if len(adult) > start_index else 0
    avg_adult_fe = sum(adult_fe[start_index:]) / len(adult_fe[start_index:]) if len(adult_fe) > start_index else 0
    avg_fertile_fe = sum(fertile_fe[start_index:]) / len(fertile_fe[start_index:]) if len(fertile_fe) > start_index else 0
    total_infection = sum(infection[start_index:]) if len(infection) > start_index else 0
    total_native_new_infections = sum(native_new_infections[start_index:]) if len(native_new_infections) > start_index else 0
    total_intro_human = sum(intro_human[start_index:]) if len(intro_human) > start_index else 0
    immune_diff = pop_immune[-1] - pop_immune[start_index] if len(pop_immune) > start_index else 0

    return (
        parsed["suppressed"],
        parsed["gen_suppressed"],
        avg_adult, avg_adult_fe, avg_fertile_fe,
        total_infection, total_intro_human, total_native_new_infections,
        parsed["mark1"], parsed["mark2"], parsed["mark3"], immune_diff
    )

def writeOutput(file, parsed, drive_number):
    headers = [
        "Generation", "Total", "Female", "Fertile Female", "Infertile Female",
        "Rate", "New Infections", "New 1", "New 2", "New 3", "New 4",
        "Pop Immune", "Intro Human", "Native New Infections"
    ]

    drive_headers = [f"{x}_{i}" for i in range(1, drive_number + 1) for x in ("WT", "DR", "R1", "R2", "R1_drive")]
    headers += drive_headers + ["R1 Introduction Year", "Female Somatic Fitness", "Conversion Rate", "ANN", "VNN"]

    df = pd.DataFrame({
        "Generation": parsed["gen"],
        "Total": parsed["total"],
        "Female": parsed["female"],
        "Fertile Female": parsed["fertile_female"],
        "Infertile Female": parsed["infertile_female"],
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
        df[f"R1_drive{i+1}"] = parsed["r1_drive"][i]
    df["Conversion Rate"] = parsed["con"]
    df["Female Somatic Fitness"] = parsed["female_fit"]
    df["R1 Rate"] = parsed["r1_rate"]
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
    clargs = "slim "
    source = args_dict.pop("source")
    for arg, val in args_dict.items():
        if isinstance(val, bool):
            clargs += f"-d {arg}={'T' if val else 'F'} "
        else:
            clargs += f"-d {arg}={val} "
    clargs += source
    return clargs.split()

def main():
    parser = ArgumentParser()
    parser.add_argument('-src', '--source', default="R1_shooter.slim", type=str,
                        help=r"SLiM file to be run.")

    parser.add_argument('-adult_females', '--NUM_ADULT_FEMALES', default=4000, type=int)
    parser.add_argument('-hv', '--H_V_RATE', default=0.4375, type=float)
    parser.add_argument('-vh', '--V_H_RATE', default=0.6, type=float)

    parser.add_argument('-drop', '--DROP', default=True)
    parser.add_argument('-drive_number', '--DRIVE_NUMBER', default=1, type=int)
    parser.add_argument('-drop_time', '--DROP_TIME', default=16, type=int)
    parser.add_argument('-drop_int', '--DROP_INTERVAL', default=52, type=int)

    parser.add_argument('-drop_1', '--DROP_1', default=0.005, type=float)
    parser.add_argument('-conversion_1', '--DRIVE_CONVERSION_RATE_1', default=0.8, type=float)
    parser.add_argument('-germline_1', '--LATE_GERMLINE_RESISTANCE_CUT_RATE_1', default=0.5, type=float)
    parser.add_argument('-embryo_1', '--EMBRYO_RESISTANCE_CUT_RATE_1', default=0.1, type=float)
    parser.add_argument('-r1_1', '--R1_OCCURRENCE_RATE_1', default=0.01, type=float)
    parser.add_argument('-m_1', '--SOMATIC_FITNESS_MULTIPLIER_M_1', default=1.0, type=float)
    parser.add_argument('-f_1', '--SOMATIC_FITNESS_MULTIPLIER_F_1', default=0.7, type=float)
    parser.add_argument('-outfile', '--output_file', default="R1_test.xlsx", type=str)
    parser.add_argument('-file', '--file1', default="Process.xlsx", type=str)

    parser.add_argument('-drop_time_r1', '--DROP_TIME_R1', default=52, type=int)
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
