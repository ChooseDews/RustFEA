import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#find the percent parallel speedup and compare to Amdahl's law
def amdahls_law(percent_parallel: float, cores: int) -> float:
    return 1 / ((1 - percent_parallel) + (percent_parallel / cores))


def get_amdahls_law_for_percent(percent_parallel: float, max_cores: int) -> list[float]:
    return [amdahls_law(percent_parallel, i) for i in range(1, max_cores + 1)]

def get_error(found_speedup: list[float], amdahls_law_speedup: list[float]) -> float:
    return (((np.array(found_speedup) - np.array(amdahls_law_speedup))) ** 2).mean()


def get_amdahls_law_for_percent_smooth(percent_parallel: float, max_cores: int):
    cores = np.linspace(1, max_cores, 100)
    return cores, [amdahls_law(percent_parallel, i) for i in cores]



def main(file='./temp/benchmark/parallel_results.csv'):
    df = pd.read_csv(file)
    one_thread_avg_time = df[df['threads'] == 1]['time'].mean()
    df['relative_time'] = df['time'] / one_thread_avg_time
    df['speedup'] = 1 / df['relative_time']


    found_speedup = df.groupby('threads')['speedup'].mean().to_list()
    print("Speedup:", found_speedup)
    max_cores = df['threads'].max()

    left, right = 0.1, 0.9
    closest_mse = float('inf')
    closest_percent = 0
    tolerance = 0.0001  # Stop when the search interval is this small
    iterations = 0
    while (right - left) > tolerance: #binary search for the percent parallel
        iterations += 1
        mid = (left + right) / 2
        mid_left = (left + mid) / 2
        mid_right = (mid + right) / 2
        mse_left = get_error(found_speedup, get_amdahls_law_for_percent(mid_left, max_cores))
        mse_mid = get_error(found_speedup, get_amdahls_law_for_percent(mid, max_cores))
        mse_right = get_error(found_speedup, get_amdahls_law_for_percent(mid_right, max_cores))
        if mse_mid < closest_mse:
            closest_mse = mse_mid
            closest_percent = mid
        if mse_left < mse_mid:
            right = mid
        elif mse_right < mse_mid:
            left = mid
        else:
            left = mid_left
            right = mid_right
    print(get_amdahls_law_for_percent(mid, max_cores))
    print(f"Closest percent: {closest_percent:.4f} with err={closest_mse:.6f} after {iterations} iterations")

    #plot the found speedup vs the amdahls law speedup for the closest percent
    amdahls_law_speedup = get_amdahls_law_for_percent(closest_percent, max_cores)
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df, x='threads', y='speedup', label='Observed Program Speedup', marker='o')

    cores, amdahls_law_speedup = get_amdahls_law_for_percent_smooth(closest_percent*1.05, int(max_cores*1.25))
    sns.lineplot(x=cores, y=amdahls_law_speedup, label=f'Amdahl\'s Law ({closest_percent * 105:.1f}%)', linestyle=':', alpha=0.5)

    cores, amdahls_law_speedup = get_amdahls_law_for_percent_smooth(closest_percent, int(max_cores*1.25))
    sns.lineplot(x=cores, y=amdahls_law_speedup, label=f'Matched Amdahl\'s Law ({closest_percent * 100:.1f}%)', linestyle='--', alpha=0.9)

    cores, amdahls_law_speedup = get_amdahls_law_for_percent_smooth(closest_percent*0.95, int(max_cores*1.25))
    sns.lineplot(x=cores, y=amdahls_law_speedup, label=f'Amdahl\'s Law ({closest_percent * 95:.1f}%)', linestyle=':', alpha=0.5)
    
    plt.title(f'Program Parallel Speedup vs Threads (MSE={closest_mse:.6f} @ {closest_percent * 100:.1f}%)')
    plt.xlabel('Threads')
    plt.ylabel('Speedup Factor')
    plt.legend()
    output = file.replace('.csv', '_plot.png')
    plt.savefig(output)
    plt.close()

if __name__ == "__main__":
    main()