import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import os
import time
import parallel_benchmark_plotter
from multiprocessing import Pool, Queue, Manager
from tqdm import tqdm
import numpy as np
import argparse


EXAMPLE_PATH = "./../examples/circle_pull_benchmark.toml"
EXAMPLE_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), EXAMPLE_PATH))
BINARY_PATH = "./../target/release/read_input"
BINARY_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), BINARY_PATH))
CMD = f"{BINARY_PATH} '{EXAMPLE_PATH}'"

DEFAULT_RESULT_ROOT = "./temp"
MIN_EXPECTED_TIME = 10 #secs
SAMPLE_SIZE = 3


n_processes = 4
RAYON_NUM_THREADS = 3

#total num cores = n_processes * RAYON_NUM_THREADS = 9


time_stamp = time.strftime("%Y-%m-%d_%H-%M-%S")
analysis_name = f"error_check_{time_stamp}"
output_dir = os.path.join(DEFAULT_RESULT_ROOT, analysis_name)
os.makedirs(output_dir, exist_ok=True)
LOG_PATH = os.path.join(output_dir, "error_count.log")
ERROR_LOG_PATH = os.path.join(output_dir, "error_log.log")
OUTPUT_FILE_PATH = os.path.join(output_dir, "error_count.csv")

def write_log(msg, log_path=LOG_PATH):
    with open(log_path, "a") as f:
        print(msg)
        f.write(f"#$ {msg}\n")

def run_cmd(cmd=CMD, env_add={}, min_time=MIN_EXPECTED_TIME):
    env_str = "; ".join([f"{key}={value}" for key, value in env_add.items()])
    write_log(f"⏳ Running command: {cmd} ⏳ | 💬 Custom environment vars: {env_str}")
    env = os.environ.copy()
    env.update(env_add)
    start_time = time.time()
    result = subprocess.run(cmd, shell=True, text=True, env=env, capture_output=True)
    dt = time.time() - start_time
    write_log(result.stdout)
    write_log(result.stderr, log_path=ERROR_LOG_PATH)
    if dt < min_time:
        raise ValueError(f"Command took {dt:.2f} seconds, expected at least {min_time} seconds. Unknown reason.")
    write_log(f"⌛ Command completed in {dt:.2f} seconds ⌛")
    failed = result.returncode != 0
    return result, dt, failed

def build_binary():
    write_log("⏳ Building binary... ⏳")
    _, dt, failed = run_cmd("cargo build --release", min_time=0.01)
    if failed:
        raise ValueError("Failed to build binary")
    write_log(f"⌛ Binary built in {dt:.2f} seconds ⌛")

def run_single_test(args):
    cmd, env = args
    _, dt, failed = run_cmd(cmd, env)
    return {
        "time": dt,
        "failed": failed
    }

def test_error_count(total_runs=400):
    os.makedirs("./temp/error_check", exist_ok=True)
    if os.path.exists(LOG_PATH): os.remove(LOG_PATH)
    if os.path.exists(ERROR_LOG_PATH): os.remove(ERROR_LOG_PATH)
    env = {
        "RUST_BACKTRACE": "1",
        "RAYON_NUM_THREADS": str(RAYON_NUM_THREADS)
    }
    build_binary()

    # Use number of CPU cores minus 1 to avoid overloading
    
    results = []
    error_count = 0
    
    # Create argument list for each process
    args_list = [(CMD, env)] * total_runs
    
    # Use tqdm with imap for progress bar
    with Pool(n_processes) as pool:
        for result in tqdm(
            pool.imap_unordered(run_single_test, args_list),
            total=total_runs,
            desc="Running tests"
        ):
            results.append(result)
            if result["failed"]:
                error_count += 1
                write_log(f"Current error count: {error_count}")
            pd.DataFrame(results).to_csv(OUTPUT_FILE_PATH, index=False)
    write_log(f"Total errors: {error_count} out of {total_runs}")

    df = pd.DataFrame(results)
    df.to_csv("./temp/error_check/error_count.csv", index=False)
    return df

def plot_error_count(df):
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    sns.countplot(x="failed", data=df, ax=axs[0])
    df = df[df["failed"] == False]
    avg_time = df["time"].mean()
    sns.histplot(x="time", data=df, ax=axs[1])
    axs[1].axvline(avg_time, color="red", linestyle="--", label=f"Average time: {avg_time:.2f} seconds")
    plt.savefig(os.path.join(output_dir, "error_count.png"))

def parse_args():
    parser = argparse.ArgumentParser(description='Run error count tests or plot existing results.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run', type=int, help='Run tests with specified number of iterations')
    group.add_argument('--plot', type=str, help='Plot results from specified CSV file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    
    if args.run:
        test_start_time = time.time()
        df = test_error_count(args.run)
        dt = time.time() - test_start_time
        print(df.describe())
        print(f"Total time: {dt:.2f} seconds")
        plot_error_count(df)
    else:  # args.plot
        df = pd.read_csv(args.plot)
        plot_error_count(df)
