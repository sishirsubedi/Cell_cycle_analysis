# import pdb


########################### simple #############################################
def match(t, p):
    count = 0
    location = {}
    # location[count]=0

    i = 0

    # pdb.set_trace()

    for j in range(len(t)):
        if i < len(p) and p[i] == t[j]:
            if count == j:
                location[count] = 0
            location[count] = location[count] + 1
            i = i + 1
        elif i == len(p):
            break
        else:
            i = 0
            count = j + 1

    return location

############################## simple ##########################################
def simplematch(t, p):
    for i in range(1, len(t)):
        j = 1
        while j <= (len(p) - 1) and t[i + j - 1] == p[j]:
            j += 1
        if j == len(p):
            print "occurs", i


############################### simple #########################################
def naive_string_search(P, T):
    m = len(T)
    n = len(P)
    indices = []
    for i in range(m - n + 1):
        if P == T[i: i + n]:
            indices.append(i)
    return indices

###################### z algorithm ############################
########################################################################


# main idea here is to find a z-box with left and right value and make an array 
# to record the maximum length of this z-box that can be generated at each index of 
# text, and if this length matches with the length of pattern then it means that 
# there is the pattern present in the text.



def combine(text,pattern):
    input = text + '$' + pattern
    return input





def zalgorithms(input):
    n = len(input)
    z = []
    for i in range(n):
        z.append(-1)
    left = right = 0

    for k in range(1, n):
        if k > right:
            left = right = k
            while right < n and input[right] == input[right - left]:
                right += 1
            z[k] = right - left
            right -= 1
        else:
            k1 = k - left
            if z[k1] < (right - k + 1):
                z[k] = z[k1]
            else:
                left = k
                while right < n and input[right] == input[right - left]:
                    right += 1
                z[k] = right - left
                right = right - 1
    return z


prime = 101

####################### Rabin Karp algorithm ############################
########################################################################

# this is the hash tagging way- instead of comparing one character you has the 
#entire pattern and then look for that hash in the text. and if that hash is present 
#then only check each characters. This way you minimize comparison.

def pattern_matching(text, pattern):
    m = len(pattern)
    n = len(text)
    pattern_hash = create_hash(pattern, m - 1)
    text_hash = create_hash(text, m - 1)

    for i in range(1, n - m + 2):
        if pattern_hash == text_hash:
            if check_equal(text[i - 1:i + m - 1], pattern[0:]) is True:
                return i - 1;
        if i < n - m + 1:
            text_hash = recalculate_hash(text, i - 1, i + m - 1, text_hash, m)
    return -1;


def check_equal(str1, str2):
    if len(str1) != len(str2):
        return False;
    i = 0
    j = 0
    for i, j in zip(str1, str2):
        if i != j:
            return False;
    return True


def create_hash(input, end):
    hash = 0
    for i in range(end + 1):
        hash = hash + ord(input[i]) * pow(prime, i)
    return hash

#Given a string representing one Unicode character, return an integer representing the Unicode
# code point of that character. For example, ord('a') returns the integer 97

def recalculate_hash(input, old_index, new_index, old_hash, pattern_len):
    new_hash = old_hash - ord(input[old_index])
    new_hash = new_hash / prime
    new_hash += ord(input[new_index]) * pow(prime, pattern_len - 1)
    return new_hash;


####################### Knuth-Morris-Pratt algorithm ############################
########################################################################

# so idea is - first you take pattern and make an array to find out if there is any repeating
#pattern inside pattern, for example in "abcxydabcef" has "abc". so now main trick is to make 
# sure when comparing this pattern to text and if there is "abc" match but "e" does not match then 
# you dont have to test "abc" in the text and you can assume "abc" is already there and start from 
# next to "abc".

def compute_temporary_array(pattern):
    n = len(pattern)
    lsp = [0 for j in range(n)]
    index = 0
    i = 1
    while i < len(pattern):
        if pattern[i] == pattern[index]:
            lsp[i] = index + 1
            index += 1
            i += 1
        else:
            if index != 0:
                index = lsp[index - 1]
            else:
                lsp[i] = 0
                i += 1
    return lsp


# KMP algorithm of pattern matching.
def kmp(text, pattern):
    lsp = compute_temporary_array(pattern)
    i = 0 # for text counter
    j = 0 #for pattern counter
    while i < len(text) and j < len(pattern):
        if text[i] == pattern[j]:
            i += 1
            j += 1
        else:
            if j != 0:
                j = lsp[j - 1]
            else:
                i += 1
    if j == len(pattern):
        return True
    else:
      return False

####################### Boyer Moore algorithm ############################
########################################################################




################################ drivers ######################################################
'''

t = "abcdabecd"
p = "abe"
str1 ="aaaab$aaaaaaaaaabaaaaaaaaa"
str2 ="abc$defghijklmnopqur"
# print match(t, p)

#simplematch(t, p)

#print naive_string_search(p, t)

#print zalgorithms(str2)


#index = pattern_matching(t,p)
#print("Index ", index)


src = 'abcxabcdabcdabcy'
sub_string = 'abcdabcy'
result = kmp(src, sub_string)
print(result)

'''
################################ run from command line #######################
#############################################################################
import sys

seq1_file = sys.argv[1]

seq2_file = sys.argv[2]

text = open(seq1_file, 'r').read()
pattern = open(seq2_file, 'r').read()

# result = kmp(text,pattern )
#print(result)


input = combine(text,pattern)
z = zalgorithms(input)
for i in z:
    if i >= len(pattern)-1:
        print 'aligns at ', i
    else:
        print ' did not find match'