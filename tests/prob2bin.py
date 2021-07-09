import sys
import argparse
import numpy as np
from scipy import sparse as sp
from lap import lapmod as g_lapmod

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n',       type=int, default=100)
    parser.add_argument('--outpath', type=str, default='data.bin')
    parser.add_argument('--seed',    type=int, default=123)
    return parser.parse_args()

def prob2bin(adj, outpath):
    n       = adj.shape[0]
    shape   = np.array(adj.shape).astype(np.int32)
    nnz     = np.array([adj.nnz]).astype(np.int32)

    indptr  = adj.indptr.astype(np.int32)
    indices = adj.indices.astype(np.int32)
    data    = adj.data.astype(np.float64)

    # _, x, y  = g_lapmod(n, data.max() - data, indptr, indices, fp_version=1)
    # mod_cost = dadj[(np.arange(n), x)].sum()
    # print(x[:10])
    # print('cost: ', mod_cost)
    
    print(f'writing {outpath}', file=sys.stderr)
    with open(outpath, 'wb') as f:
        _ = f.write(bytearray(shape))
        _ = f.write(bytearray(nnz))
        _ = f.write(bytearray(indptr))
        _ = f.write(bytearray(indices))
        _ = f.write(bytearray(data))
        f.flush()

if __name__ == "__main__":
    args = parse_args()
    np.random.seed(args.seed)

    dadj    = np.random.uniform(0, 1, (args.n, args.n))
    m       = np.random.uniform(0, 1, dadj.shape) > 0.5
    dadj[m] = 0
    adj     = sp.csr_matrix(dadj)

    prob2bin(adj, args.outpath)
    
    print('-' * 50)