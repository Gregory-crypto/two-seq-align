

# supplementary functions

def isnt_comp(primer, max_comp):

    comp = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    rev_part = ''

    for i in primer[-1:-1 - max_comp:-1]:
        rev_part += comp[i]

    for i in range(len(primer) - max_comp + 1):
        if rev_part == primer[i:i + max_comp]:
            return False
        

    return True

def gc_check(content, min_gc, max_gc, prim_leng):
    if min_gc <= ((content['G'] + content['C']) / prim_leng) * 100 <= max_gc:
        return True
    else:
        return False

def tm_check(content, min_t, max_t, prim_leng):
    if min_t <= (64.9 + 41 * (content['G'] + content['C'] - 16.4) / prim_leng) <= max_t:
        return True
    else:
        return False

    

def get_primer(seq, prim_leng, min_gc = 50, max_gc = 60, min_t = 50, max_t = 60, max_comp = 4):
    
    # preprocessing
    seq = seq.upper()
    seq_length = len(seq)

    # first step

    content = {
        'A': 0,
        'T': 0,
        'G': 0,
        'C': 0
    }

    for i in seq[0:prim_leng]:
        content[i] += 1
    
    gc_content = gc_check(content, min_gc, max_gc, prim_leng)
    tm = tm_check(content, min_t, max_t, prim_leng)
    is_not_comp = isnt_comp(seq[0:prim_leng], max_comp + 1)


    if all([gc_content, tm, is_not_comp]):
        return seq[0:prim_leng]


    for i in range(seq_length - prim_leng):
        # content editing
        content[seq[i]] -= 1
        content[seq[prim_leng + i]] += 1

        gc_content = gc_check(content, min_gc, max_gc, prim_leng)
        tm = tm_check(content, min_t, max_t, prim_leng)
        is_not_comp = isnt_comp(seq[i + 1:prim_leng + i + 1], max_comp + 1)


        if all([gc_content, tm, is_not_comp]):
            return seq[i + 1:prim_leng + i + 1]
        
    return False

# returns 'AGCATGCATCGATCGACTAGC'
# seq = 'TAGCATGCATCGATCGACTAGCTACGATCGATCGACTAATTACTACGGCCGCGATCGACCGTACTAATCGATCATGTAATATTACGATCGAT'
# prim_leng = 21

# returns False
# seq = 'ATATATATCGATGCTATATGCGCTATATACTGACTAGCATCGATCGATATAAAA'
# prim_leng = 20

# print(get_primer(seq, prim_leng))
