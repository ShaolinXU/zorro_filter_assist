#!/Users/xushaolin/miniconda3/bin/python

from Bio import AlignIO
import numpy
import sys
import getopt

arguments = sys.argv[1:]

try:  
    opts, args = getopt.getopt(arguments, "-s:-m:-o:-t:-h",["sequence","mask","output","threshold","help"])
except getopt.GetoptError:  
    # print help information and exit: 
    print("Your options have some error")
    print('''Modified by Shaolin.Xu (2020-07-01) from script from Schwentner, Martin (DOI:https://doi.org/10.5061/)
        -s for input sequence, require string
        -m for mask file, require string
        -t for threshold, require string
        -o for output file, require string
        -h (--help) for help info''')
    exit()

for opt, arg in opts:
    if opt in ['-h','--help']:
        print('''Modified by Shaolin.Xu (2020-07-01) from script from Schwentner, Martin (DOI:https://doi.org/10.5061/)
        -s (--sequence) for input sequence, require string
        -m (--mask) for mask file, require string
        -t (--threshold) for threshold, require string
        -o (--output) for output file, require string
        -h (--help) for help info''')
        exit()
    elif opt in ['-s']:
        seq = arg
    elif opt in ['-m']:
        mask = arg
    elif opt in ['-o']:
        out = arg
    elif opt in ['-t']:
        threshold = float(arg)

#threshold = 5 #float(raw_input('What HMM value do you which to exclude? '))

alignment = AlignIO.read(seq, 'fasta')
align_array = numpy.array([list(rec) for rec in alignment], numpy.character, order="F")
zorro_handle = open(mask, 'r')
count = 0
cut = align_array[:,:1]
for j in zorro_handle:
        if float(j.rstrip('\n')) >= threshold:
                cut = numpy.column_stack((cut, align_array[:,count]))
        count += 1
cut = cut[:,1:]
zorro_handle.close()

outfile = open(out, 'w')
for q in range(len(alignment)):
        seq = ''
        for p in cut[q,:]:
                seq += bytes.decode(p)
        outfile.write('>' + alignment[q].id + '\n' + seq + '\n')
outfile.close()
