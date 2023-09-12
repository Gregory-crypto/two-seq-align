import re

def find_target(seq, pam, pam_ori, len_tar):

    # storage for targets

    target_list = []

    # processing pam

    reg_ipam = ''
    reg_iipam = ''

    base_dict = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'R': '(G|A)',
        'Y': '(T|C)',
        'K': '(G|T)',
        'M': '(A|C)',
        'S': '(G|C)',
        'W': '(A|T)',
        'B': '(G|T|C)',
        'D': '(G|A|T)',
        'H': '(A|C|T)',
        'V': '(G|C|A)',
        'N': '(A|G|C|T)',

    }

    for i in pam:
        reg_ipam += base_dict[i]

    
    # complementary string

    comp = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }

    rev_seq = ''

    for i in seq[::-1]:
        rev_seq += comp[i]

    # find pam

    # if pam orientation is R

    if pam_ori == 'R':

        while True:
            match = re.search(reg_ipam, seq)
            if match:
                if len(seq[:match.start() - 1]) >= len_tar:
                    target_list.append(seq[match.start() - len_tar:match.start()])
            else:
                break
            seq = seq[match.start() + 1:]

        while True:
            match = re.search(reg_ipam, rev_seq)
            if match:
                if len(rev_seq[:match.start() - 1]) >= len_tar:
                    target_list.append(rev_seq[match.start() - len_tar:match.start()])
            else:
                break
            rev_seq = rev_seq[match.start() + 1:]
        
    # if pam orientation is L

    else:

        while True:
            match = re.search(reg_ipam, seq)
            if match:
                if len(seq[match.end():]) >= len_tar:
                    target_list.append(seq[match.end():match.end() + len_tar])
            else:
                break
            seq = seq[match.start() + 1:]

        while True:
            match = re.search(reg_ipam, rev_seq)
            if match:
                if len(rev_seq[match.end():]) >= len_tar:
                    target_list.append(rev_seq[match.end():match.end() + len_tar])
            else:
                break
            rev_seq = rev_seq[match.start() + 1:]


    return target_list


# seq = "TAGCTACGATCGATCGTTTCTAGCTACGATGCAAGAAAGATCGATCGATCGACGTACG"
# pam = "YTTN"
# pam_ori = "L"
# target_length = 21

# print(find_target(seq, pam, pam_ori, target_length))