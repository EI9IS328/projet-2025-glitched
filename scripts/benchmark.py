from io import TextIOWrapper
import re
import subprocess
import shutil
import sys
import os
import itertools


NUMBER_REGEX_PART = r"(\d+(?:.\d+)?(?:e[+-]?\d+)?)"
TIME_REGEX_PART = rf"{NUMBER_REGEX_PART} seconds\."
BYTES_REGEX_PART = rf"{NUMBER_REGEX_PART} Bytes\."
TIME_TOTAL_REGEX = rf"Time Kernel Total : {TIME_REGEX_PART}"
TIME_KERNEL_REGEX = rf"Time Spent Simulating : {TIME_REGEX_PART}"
TIME_SNAPSHOTS_REGEX = rf"Time Making and Saving Snapshots : {TIME_REGEX_PART}"
TIME_SISMOS_REGEX = rf"Time Making Sismos : {TIME_REGEX_PART}"
TIME_OUTPUT_SISMOS_REGEX = rf"Time Saving Outputs : {TIME_REGEX_PART}"
TOTAL_BYTES_REGEX = rf"Total written data: {BYTES_REGEX_PART}"
BYTES_RECEIVERS_REGEX = rf"receiver_\d+: {BYTES_REGEX_PART}"
SNAPSHOT_RECEIVERS_REGEX = rf"snapshot\d+.bin: {BYTES_REGEX_PART}"

print(TIME_TOTAL_REGEX)


def benchmark(
    bin: str,
    output_dir: str,
    output_file: TextIOWrapper,
    exs: list[int],
    eys: list[int],
    ezs: list[int],
    orders: list[int],
    formats_sismos: list[str],
    formats_snaps: list[str],
    receivers: list[list[tuple[int, int, int]]],
):
    output_file.write(
        "ex|ey|ez|order|format_sismos|format_snaps|receivers|global|kernel|make_snapshots|make_sismos|output_sismos|total_bytes|total_snap_bytes|total_sismos_bytes\n"
    )
    tmp_path = os.path.join(output_dir, "tmp")
    if os.path.exists(tmp_path):
        shutil.rmtree(tmp_path)
    os.makedirs(tmp_path, exist_ok=True)
    combs = [exs, eys, ezs, orders, formats_sismos, formats_snaps, receivers]
    combs_product = list(itertools.product(*combs))
    done = 0
    for comb in combs_product:
        # we set up everything
        ex, ey, ez, order, format_sismos, format_snaps, cur_receivers = comb
        param_args = f"--ex {ex} --ey {ey} --ez {ez} -o {order}"
        static_args = "--implem makutu"
        output_receivers = os.path.join(output_dir, "tmp", f"receiver_{done}")
        receivers_file_path = os.path.join(output_dir, "tmp", f"receivers_{done}.txt")
        snaps = os.path.join(output_dir, "tmp", f"snaps_{done}")
        measures = os.path.join(output_dir, f"measures_{done}")
        files_args = f"--watched-receivers {receivers_file_path} --output-receivers {output_receivers} --output-receivers-format {format_sismos} --snapshot-folder-path {snaps} --output-measures {measures} --snapshot-format {format_snaps}"
        cmdline = f"{bin} {param_args} {static_args} {files_args}"
        os.mkdir(snaps)

        # we create a temporary file for handling receivers
        with open(receivers_file_path, "w+") as receivers_file:
            for receiver in cur_receivers:
                receivers_file.write(f"{receiver[0]};{receiver[1]};{receiver[2]}\n")
        print(f"running exp {done + 1}/{len(combs_product)} ({cmdline})")
        process_rez = subprocess.run(
            cmdline, shell=True, capture_output=True, text=True
        )
        if process_rez.returncode != 0:
            shutil.rmtree(tmp_path)
            print(process_rez.stderr, file=sys.stderr)
            exit(1)
        os.remove(receivers_file_path)
        os.remove(output_receivers)
        shutil.rmtree(snaps)

        # TODO: make this safer
        with open(measures, "r") as measures_file:
            measures_file_content = measures_file.read()
            try:
                time_total = re.findall(TIME_TOTAL_REGEX, measures_file_content)[0]
                time_kernel = re.findall(TIME_KERNEL_REGEX, measures_file_content)[0]
                time_snaps = re.findall(TIME_SNAPSHOTS_REGEX, measures_file_content)[0]
                time_making_sismos = re.findall(
                    TIME_SISMOS_REGEX, measures_file_content
                )[0]
                time_saving_sismos = re.findall(
                    TIME_OUTPUT_SISMOS_REGEX, measures_file_content
                )[0]
                total_written_data = re.findall(
                    TOTAL_BYTES_REGEX, measures_file_content
                )[0]
                written_receiver_data = re.findall(
                    BYTES_RECEIVERS_REGEX, measures_file_content
                )[0]
                snapshots_receiver_data = re.findall(
                    SNAPSHOT_RECEIVERS_REGEX, measures_file_content
                )
                total_snap_bytes = sum(map(lambda x: int(x), snapshots_receiver_data))
                csv_line = f"{ex}|{ey}|{ez}|{order}|{format_sismos}|{format_snaps}|{len(receivers[0])}|{time_total}|{time_kernel}|{time_snaps}|{time_making_sismos}|{time_saving_sismos}|{total_written_data}|{total_snap_bytes}|{written_receiver_data}\n"
                output_file.write(csv_line)
            except Exception as e:
                print(f"Exception occured: {e}", file=sys.stderr)
                print(measures_file_content)
                shutil.rmtree(tmp_path)
                break
        os.remove(measures)
        done += 1


if __name__ == "__main__":
    OUTPUT_DIR = "./exp1"
    OUTPUT_FILE = "result.csv"
    EXS = [10]
    EYS = [5]
    EZS = [5]
    ORDERS = [1]
    FORMATS_SISMOS = ["bin"]
    FORMATS_SNAPS = ["bin", "plain"]
    RECEIVERS = [
        [(1, 4, 3), (122, 312, 451)],
    ]
    BIN = "./build/debug/bin/semproxy"

    COMPILE = True

    COMPILE_COMMAND = "cmake --build ./build/debug"

    if COMPILE:
        compile_rez = subprocess.run(COMPILE_COMMAND, shell=True, capture_output=True)
        if compile_rez.returncode == 1:
            print(
                f"Compilation failed:\n\n{compile_rez.stderr.decode('utf-8')}",
                file=sys.stderr,
            )
            exit(1)

    output_file_path = os.path.join(OUTPUT_DIR, OUTPUT_FILE)
    with open(output_file_path, "w+") as output_file:
        benchmark(
            BIN,
            OUTPUT_DIR,
            output_file,
            EXS,
            EYS,
            EZS,
            ORDERS,
            FORMATS_SISMOS,
            FORMATS_SNAPS,
            RECEIVERS,
        )
    print(f"done. results in {output_file_path}")
