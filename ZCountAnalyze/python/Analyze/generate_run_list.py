# generate_run_list.py

def extract_run_numbers(csv_path: str, output_path: str):
    run_numbers = set()

    with open(csv_path, "r") as f:
        for line in f:
            if ":" in line:
                run = line.split(":")[0].strip()
                if run.isdigit():
                    run_numbers.add(run)

    with open(output_path, "w") as out:
        for run in sorted(run_numbers, key=int):
            out.write(run + "\n")

# Usage
extract_run_numbers("/eos/cms/store/group/comm_luminosity/ZCounting/2022/brilcalcByLS/BRIL/Golden/byLS_2022_355100_362760_Golden.csv", "Golden_run_list.txt")
