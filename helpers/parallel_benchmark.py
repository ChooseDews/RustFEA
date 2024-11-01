import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import os
import time
import parallel_benchmark_plotter


EXAMPLE_PATH = "./../examples/circle_pull_benchmark.toml"
EXAMPLE_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), EXAMPLE_PATH))
BINARY_PATH = "./../target/release/read_input"
BINARY_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), BINARY_PATH))
CMD = f"{BINARY_PATH} '{EXAMPLE_PATH}'"


timestamp = time.time()
folder_name = f"parallel_benchmark_{timestamp}"
OUTPUT_DIR = f"./temp/benchmark/"
OUTPUT_DIR = os.path.join(OUTPUT_DIR, folder_name)
os.makedirs(OUTPUT_DIR, exist_ok=True)
LOG_PATH = f"{OUTPUT_DIR}/parallel_benchmark.log"
MIN_EXPECTED_TIME = 10 #secs
SAMPLE_SIZE = 1



def write_log(msg):
    with open(LOG_PATH, "a") as f:
        print(msg)
        f.write(f"#$ {msg}\n")

def run_cmd(cmd=CMD, env_add={}, min_time=MIN_EXPECTED_TIME):
    env_str = "; ".join([f"{key}={value}" for key, value in env_add.items()])
    write_log(f"⏳ Running command: {cmd} ⏳ | 💬 Custom environment vars: {env_str}")
    env = os.environ.copy()
    env.update(env_add)
    start_time = time.time()
    result = subprocess.run(cmd, shell=True, check=True, text=True, env=env)
    dt = time.time() - start_time
    if dt < min_time:
        raise ValueError(f"Command took {dt:.2f} seconds, expected at least {min_time} seconds. Unknown reason.")
    write_log(f"⌛ Command completed in {dt:.2f} seconds ⌛")
    return result, dt

def build_binary():
    write_log("⏳ Building binary... ⏳")
    _, dt = run_cmd("cargo build --release", min_time=0.01)
    write_log(f"⌛ Binary built in {dt:.2f} seconds ⌛")

def test_time_over_thread_count(to_thread=10):
    os.makedirs("./temp/benchmark", exist_ok=True)
    if os.path.exists(LOG_PATH): os.remove(LOG_PATH)
    env = {
        "RUST_BACKTRACE": "full"
    }
    build_binary()
    thread_counts = list(range(1, min(to_thread, os.cpu_count()) + 1))
    print(f"Testing with thread counts: {thread_counts}")
    results = []
    
    for i, threads in enumerate(thread_counts):
        env["RAYON_NUM_THREADS"] = str(threads)
        print(f"Running with {threads} threads")
        
        for _ in range(SAMPLE_SIZE):  # 3 runs per thread count
            for attempt in range(3):
                try:
                    _, dt = run_cmd(CMD, env)
                    results.append({
                    'threads': threads,
                        'time': dt,
                    })
                    break
                except Exception as e:
                    write_log(f"Error: {e}")
                    if attempt == 2:
                        write_log(f"Failed after 3 attempts. Raising error.")
                        raise e
        if i == 0:
            continue


        df = pd.DataFrame(results)
        output_path = os.path.join(OUTPUT_DIR, "parallel_results.csv")
        df.to_csv(output_path, index=False)
        parallel_benchmark_plotter.main(output_path)



    return df


if __name__ == "__main__":
    test_start_time = time.time()
    df = test_time_over_thread_count()
    dt = time.time() - test_start_time
    print(df)
    print(f"Total time: {dt:.2f} seconds")
