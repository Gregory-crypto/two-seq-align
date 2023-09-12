## `Align()` function

The code is a Python function called `align()`, which performs **Needleman-Wunsch sequence alignment** on two input sequences. Needleman-Wunsch is a global alignment algorithm, which means that it finds the optimal alignment between two sequences over their entire lengths.

The `align()` function takes two sequences as input, `seqi` and `seqii`, and returns a list of two aligned sequences, `al_seqi` and `al_seqii`. The aligned sequences are returned as strings, with any gaps indicated by the `-` character.

The `align()` function works by first creating a scoring matrix, `M`. The scoring matrix is a table that stores the optimal alignment score for each possible alignment of the two input sequences. The scores in the matrix are calculated using the following parameters:

- `match`: The score for aligning two matching residues
- `mismatch`: The score for aligning two mismatched residues
- `open_gap`: The penalty for opening a new gap in the alignment
- `gap`: The penalty for extending an existing gap in the alignment

The scoring matrix is initialized with the following values:
- `M[0, 1] += open_gap`
- `M[1, 0] += open_gap`
- `M[1:, 0] = np.linspace(-2, -2 + (ni - 1) * gap, ni)`
- `M[0, 1:] = np.linspace(-2, -2 + (nii - 1) * gap, nii)`

The first two values indicate the penalty for opening a new gap at the beginning of either sequence. The third and fourth values indicate the penalty for extending a gap that is already present in either sequence.

Once the scoring matrix has been initialized, the `align()` function iterates over the matrix, filling in the remaining values. The value for each cell in the matrix is calculated as the maximum of the following three values:
- The score for aligning the current residues in the two input sequences, plus the score for the optimal alignment of the remaining residues
- The score for extending a gap in the first input sequence, plus the score for the optimal alignment of the remaining residues in the second input sequence
- The score for extending a gap in the second input sequence, plus the score for the optimal alignment of the remaining residues in the first input sequence

Once the scoring matrix has been filled in, the `align()` function traces back from the bottom right corner of the matrix to find the optimal alignment. The tracing process starts by adding the last residue of each input sequence to the aligned sequences. Then, the algorithm recursively traces back through the scoring matrix, adding the next residue from each input sequence to the aligned sequences, depending on the value of the scoring matrix cell.

The `align()` function terminates when it reaches the top left corner of the scoring matrix. At that point, the aligned sequences are complete.

Here is an example of how to use the `align()` function:

```py
import align

seq1 = 'CACACAGTGACTAGCTAGCTACGATC'
seq2 = 'CACACAGTCGACTAGCTAGCACGATC'

aligned_sequences = align(seq1, seq2)

print(aligned_sequences)
```

```py
Output:

['CACACAGTGACTAGCTAGCTACGATC', 'CACACAGTGACTAGCTAGCACGATC']
```
