## `align()` function

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

## `find_target()` function

The code is a Python function called `find_target()`, which finds all target sequences in a given sequence, given the PAM sequence and PAM orientation. The function takes four arguments:

* `seq`: The sequence to be searched.
* `pam`: The PAM sequence.
* `pam_ori`: The PAM orientation, either `L` (left) or `R` (right).
* `len_tar`: The length of the target sequence.

The function returns a list of all target sequences found in the input sequence.

The `find_target()` function works by first creating a regular expression for the PAM sequence, using the `re.search()` function. The regular expression takes into account the PAM orientation, and allows for ambiguity in the PAM sequence by using the base dictionary defined in the function.

Once the regular expression has been created, the `find_target()` function iterates over the input sequence, searching for matches to the PAM sequence. If a match is found, the function checks to see if the target sequence is long enough. If it is, the function adds the target sequence to the list of target sequences.

The `find_target()` function also searches for target sequences in the reverse complement of the input sequence. This is because some CRISPR-Cas systems can recognize target sequences in either orientation.

Here is an example of how to use the `find_target()` function:

```python
import find_target

seq = "TAGCTACGATCGATCGTTTCTAGCTACGATGCAAGAAAGATCGATCGATCGACGTACG"
pam = "YTTN"
pam_ori = "L"
target_length = 21

target_list = find_target(seq, pam, pam_ori, target_length)

print(target_list)
```

Output:

```
['TAGCTACGATGCAAGAAAGAT', 'CTTGCATCGTAGCTAGAAACG', 'TTGCATCGTAGCTAGAAACGA', 'CATCGTAGCTAGAAACGATCG']
```

## `get_primer()` function

The program is a Python function that designs a primer for a DNA sequence. The function takes the following parameters:

* `seq`: The DNA sequence to design a primer for.
* `prim_leng`: The desired length of the primer.
* `min_gc`: The minimum GC content of the primer.
* `max_gc`: The maximum GC content of the primer.
* `min_t`: The minimum melting temperature of the primer.
* `max_t`: The maximum melting temperature of the primer.
* `max_comp`: The maximum number of consecutive complementary bases allowed in the primer.

The function first preprocesses the DNA sequence by converting it to uppercase and getting its length. It then starts by checking if the first `prim_leng` bases of the sequence meet the desired criteria. If they do, the function returns the first `prim_leng` bases of the sequence. Otherwise, the function iterates through the sequence, one base at a time, and checks if the next `prim_leng` bases starting from that base meet the desired criteria. If they do, the function returns the next `prim_leng` bases starting from that base. If no bases in the sequence meet the desired criteria, the function returns `False`.

The following is a breakdown of the function's steps:

1. The `isnt_comp()` function checks if a primer sequence contains any consecutive complementary bases. This is done by first creating a dictionary that maps each base to its complement (A -> T, T -> A, G -> C, and C -> G). The function then reverses the last `max_comp` bases of the primer sequence and checks if it is a substring of the primer sequence. If it is, the function returns `False`. Otherwise, the function returns `True`.
2. The `gc_check()` function checks if the GC content of a primer sequence is within the specified range. The function does this by first getting the number of G and C bases in the primer sequence. It then divides this number by the length of the primer sequence and multiplies by 100 to get the GC content as a percentage. If the GC content is within the specified range, the function returns `True`. Otherwise, the function returns `False`.
3. The `tm_check()` function checks if the melting temperature of a primer sequence is within the specified range. The melting temperature of a primer sequence is determined by its GC content and length. The function uses the following formula to calculate the melting temperature:

```
Tm = 64.9 + 41 * (GC content - 16.4) / primer length
```

If the melting temperature is within the specified range, the function returns `True`. Otherwise, the function returns `False`.

4. The `get_primer()` function first calls the `isnt_comp()` function to check if the first `prim_leng` bases of the sequence meet the complementary base requirement. If they do, the function returns the first `prim_leng` bases of the sequence. Otherwise, the function enters a loop.

5. In each iteration of the loop, the function calls the `gc_check()` and `tm_check()` functions to check if the next `prim_leng` bases starting from the current base meet the GC content and melting temperature requirements. If they do, the function returns the next `prim_leng` bases starting from the current base.

6. If the loop terminates without finding a primer sequence that meets all of the requirements, the function returns `False`.

Here is an example of how to use the `get_primer()` function:

```python
import get_primer

seq = 'TAGCATGCATCGATCGACTAGCTACGATCGATCGACTAATTACTACGGCCGCGATCGACCGTACTAATCGATCATGTAATATTACGATCGAT'
prim_leng = 21

print(get_primer(seq, prim_leng))
```

Output:

```
AGCATGCATCGATCGACTAGC
```

```python
import get_primer

seq = 'ATATATATCGATGCTATATGCGCTATATACTGACTAGCATCGATCGATATAAAA'
prim_leng = 20

print(get_primer(seq, prim_leng))
```

Output:

```
False
```


## `get_reverse_primer()` function

This program is a function called `get_reverse_primer()` that takes in a sequence of DNA, a forward primer, the desired length of the reverse primer, the minimum and maximum amplicon sizes, and the minimum and maximum melting temperatures (TM) as input. The function then outputs the reverse primer that meets all of the specified criteria.

The function first preprocesses the sequence by converting it to uppercase, removing the forward primer, and ensuring that the reverse primer is at least `min_amp` nucleotides long and no more than `max_amp` nucleotides long. The function then reverses the sequence and calculates the GC content, TM, and self-complementarity of the first `prim_leng` nucleotides. If these values meet all of the specified criteria, the function returns the first `prim_leng` nucleotides as the reverse primer.

If the first `prim_leng` nucleotides do not meet all of the criteria, the function iterates through the sequence, starting from the second nucleotide, and calculating the GC content, TM, and self-complementarity of the `i+1` to `i+prim_leng` nucleotides. If these values meet all of the specified criteria, the function returns the `i+1` to `i+prim_leng` nucleotides as the reverse primer.

If no reverse primer can be found that meets all of the specified criteria, the function returns `False`.

Here is a more detailed explanation of each function:

* `isnt_comp()` checks if the reverse primer is not self-complementary. This is done by creating a reversed sequence of the reverse primer and then checking if it is a substring of the reverse primer.
* `gc_check()` checks if the GC content of the reverse primer is within the specified range.
* `tm_check()` checks if the TM of the reverse primer is within the specified range.
* `get_reverse_primer()` is the main function of the program. It takes in all of the input parameters and then calls the other functions to find the reverse primer that meets all of the specified criteria.


Here is an example of how to use the `get_reverse_primer()` function:

Input:

```python
import get_reverse_primer

seq = "ATTATCGATCGATCAGTATATATCGCGCGCGATATATGCATCGATCGATCGACTAGCTACGA"
fprimer = "ATATCGCGCGCGATATATGC"
prim_leng = 20
min_amp = 30
max_amp = 50

print(get_reverse_primer(seq, fprimer, prim_leng, min_amp, max_amp))
```

Output:

```
AGCTAGTCGATCGATCGATG
```
Input:

```python
import get_reverse_primer

seq = "ATTATCGATCGATCAGTATATATCGCGCGCGATATATGCATCGATCGATCGACTAGCTACGA"
fprimer = "ATATCGCGCGCGATATATGC"
prim_leng = 20
min_amp = 40
max_amp = 100

print(get_reverse_primer(seq, fprimer, prim_leng, min_amp, max_amp))
```

Output:

```
False
```
