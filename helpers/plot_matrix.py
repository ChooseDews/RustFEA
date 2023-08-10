import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps
import argparse

def plot_matrix(file_name, output_image):
    with open(file_name,"r") as file:
        row = []
        col = []
        value = []
        for line in file:
            line = line.strip()
            line = line.split(",")
            row.append(int(line[0]))
            col.append(int(line[1]))
            value.append(float(line[2]))

    max_row = max(row)
    max_col = max(col)
    assert max_row == max_col
    print(f"max_row = {max_row}")
    K = sps.coo_matrix((value, (row, col)), shape=(max_row+1, max_row+1)).toarray()
    plt.spy(K)
    plt.savefig(output_image)
    #check if K is symmetric
    tol = 1e-5
    assert np.allclose(K, K.T, atol=tol)
    #check if K is positive definite within a tolerance
    tol = 1e-3
    w, v = np.linalg.eig(K)
    print(f"eigenvalues = {w}")
    #print all negative eigenvalues
    print(f"negative eigenvalues = {w[w < 0]}")
    assert np.all(w > -tol)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Load matrix from txt and plot it.")
    parser.add_argument("--file", type=str, default="matrix.txt", help="Name of the matrix file to be loaded.")
    parser.add_argument("--output", type=str, default="matrix.png", help="Name of the output image file.")

    args = parser.parse_args()

    plot_matrix(args.file, args.output)