#some operations on a dictionary
#definition of a dictionary
telephone = {'rob': 4098, 'lidia': 4139,'hank':4325}
#adding a new entry
telephone['harry'] = 4100
#changing an existing one
telephone['lidia'] = 4127
#deleting an entry
del telephone['hank']
#show a textual representation of the dictionary
print(telephone)

#>>> {'rob': 4098, 'harry': 4100, 'lidia': 4127}
#which are the keys in the dictionary?
print(telephone.keys())

#>>> dict_keys(['rob', 'harry', 'lidia'])

#which are the values in the dictionary?
print(telephone.values())

# is the entry a key in the dictionary?
print('lidia' in telephone)


#! /usr/bin/python


def readTextFile(file):
    handle = open(file)
    lines = []
    for line in handle:
        line=line.strip()
        lines.append(line)
    handle.close()
    return lines

def readFasta(file):
    lines=readTextFile(file)
    lines=lines[1:]
    result = ''.join(lines)
    return result.lower()

def retrieveCodons(seq,start):
    result = []
    stop = ['taa','tga','tag']

    while (start+3)<=len(seq):
        codon=seq[start:(start+3)]
        result.append(codon);
        if(codon in stop):
            start=len(seq);
        else:
            start+=3;
    return result;

def translate(map,seq):
    result = ''
    start = 0;
    for codon in retrieveCodons(seq,0):
        result+=map[codon];
    return result

file = 'codonCodingTable.txt'
with open(file) as handle:
	map={}
	for line in handle:
    		map[line[:3]]=line[4]

sequence= translate(map,readFasta('sequence.fasta'))
with open("translation.txt", 'w') as handle:
	handle.write(sequence)

import re
pattern="([\w.-]+)@([\w.-]+)"
string="You can contact Alice through alice-b@google.com about her dishwasher"
match =re.search(pattern, string)
print(match.group())
print(match.group(1))
print(match.group(2))

import re
def testRE(pattern,string):
	match=re.search(pattern,string)
	if match!=None:
		print(string + " matches {0}".format(match.group(1)))
	else:
		print(string + " matches {0}".format(match))

strings=[]
strings.append("1234")
strings.append("-1234")
strings.append("this is an integer -1234.")
strings.append("Is this an integer -1234?")
strings.append("this is a float -1.234.")
strings.append("this is a nice 1.4 float.")
strings.append("this is an important integer 10 to catch.")
strings.append("this is not an integer -10Bis.")
pattern=r"(?:^|\s)(-?\d+)[.?!]?(?:\s|$)"
#(?: leaves the grouping by the parenthesis, but removes the capturing
for string in strings:
	testRE(pattern,string)

    dict={}
bases=list("actg")
for first in bases:
	for second in bases:
		for third in bases:
			dict[first+second+third]=0
print(dict)
