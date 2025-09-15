# Dense Matrix–Vector Multiply (Sequential, 1‑D MPI, 2‑D MPI)

This package contains three reference implementations and notes that match your assignment:

- `seq_matvec.c` – Sequential baseline.
- `mpi_matvec_1d.c` – Row‑wise 1‑D block partitioning using MPI.
- `mpi_matvec_2d.c` – 2‑D block partitioning on a `q × q` process grid using MPI (requires `p = q^2` and `n % q == 0`).

All codes initialize `A[i,j] = i + j + 1` and `x[j] = 1`. They print the first 10 entries of `y = A x` for quick sanity checking.

## Build

```bash
# Sequential
gcc -O2 seq_matvec.c -o seq_mv.exe

# MPI 1-D
mpicc -O2 mpi_matvec_1d.c -o mpi_mv1d

# MPI 2-D
mpicc -O2 mpi_matvec_2d.c -o mpi_mv2d -lm
```

## Run

```bash
# Sequential (n = 1024 for example)
./seq_mv.exe 1024

# 1-D rowwise with p processes (e.g., p=4)
mpirun -np 4 ./mpi_mv1d.exe 1024

# 2-D block with p = q*q (e.g., p=9 -> q=3) and n divisible by q
mpicc -O2 mpi_matvec_2d.c -o mpi_mv2d -lm
mpirun -np 9 ./mpi_mv2d 1026

mpirun -np 9 ./mpi_mv2d.exe 1023   # will error (since 1023 % 3 != 0)
mpirun -np 9 ./mpi_mv2d.exe 1023   # choose n divisible by 3 instead, e.g., 1026
mpirun -np 9 ./mpi_mv2d.exe 1026
```

## Notes

- The 1‑D code uses `MPI_Bcast` for the vector `x`, `MPI_Scatterv` for rows of `A`, and `MPI_Gatherv` for the pieces of `y`.
- The 2‑D code forms a Cartesian grid, distributes `A` as contiguous `b × b` blocks via subarray datatypes and point‑to‑point sends, scatters `x` blocks along the top row, broadcasts them down each column, computes local partials, reduces across columns, then gathers the `y` blocks from column‑0 processes.
- You can adapt the same structure to other languages (C++, Fortran, Python/mpi4py).

## α–β–γ Cost Model (summary)

Let `α` be latency (per message), `β` the time per double sent, and `γ` time per flop.

- **Sequential:** `T_seq ≈ 2γ n^2` flops.

- **1‑D (rowwise):**  
  Communication: broadcast `x` (`α log p + β n`) + gather `y` (`α log p + β n`).  
  Compute: `2γ n^2 / p`.  
  **Total:** `T_1D ≈ 2γ n^2 / p + 2α log p + 2β n` (plus small constants).

- **2‑D (q × q grid, p=q^2):**  
  Communication: `x` broadcast down columns (`α log q + β n/q`) + reduction of partial `y` along rows (`α log q + β n/q`) + (optional) final gather of `y` to root (`α log p + β n`).  
  Compute: `2γ n^2 / p`.  
  **Total (distributed y):** `T_2D ≈ 2γ n^2 / p + 2α log q + 2β (n/q)`;  
  **Total (gathered y):** add `α log p + β n` if collecting to a single root.

For large `p`, the 2‑D scheme reduces bandwidth from `Θ(n)` (1‑D) to `Θ(n/√p)`, giving superior scalability.