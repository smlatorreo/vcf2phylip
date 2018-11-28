#!/usr/bin/env python3
__author__ = "Sergio Latorre"
__license__ = "GPLv3.0"
__version__ = "1.0"
__email__ = "sergio.latorre at tuebingen.mpg.de"

from sys import argv

if '-i' not in argv or '-h' in argv:
    print('\nUSAGE: vcf2phylip.py <options> -i <input.vcf>')
    print('\n-f: Apply positions tagged with anything but PASS')
    print('-m <float>: Minimum allele frequency threshold expressed as fraction')
    print('-h: Prints this message')
    print('\nSergio Latorre\nsergio.latorre \'at\' tuebingen.mpg.de\n')
    exit()

if '-f' in argv:
    filt = True
else:
    filt = False

inp = argv[argv.index('-i') + 1]
# Create list with samples
headerlines = 0
with open(inp, 'r') as vcf:
    for line in vcf.readlines():
        if line.startswith('#CHROM'):
           samples = line.rsplit()[9:]
           headerlines += 1
        elif line.startswith('#'):
            headerlines += 1

for sample in samples:
    if len(sample) > 10:
        print('Please make sure the sequence {} is up to 10 characters'.format(sample))
        exit()

nsamples = len(samples)

if '-m' in argv:
    maf = float(argv[argv.index('-m') + 1])
else:
    maf = 0
maf = int(nsamples * maf)

# Main function
def parser(line, nsamples, maf):
    REF = line.rsplit()[3]
    ALT = line.rsplit()[4]
    lineout = ''
    for sample in range(9, 9+nsamples):
        if line.rsplit()[sample][0] == '0':
            lineout += REF
        elif line.rsplit()[sample][0] == '1':
            lineout += ALT
        else:
            lineout += 'N'
    if min((lineout.count(REF)), (lineout.count(ALT))) <= maf:
        lineout = ''
    return lineout

# Checking filter flag and passing lines to the main function
output = {sample:'' for sample in range(nsamples)}
with open(inp, 'r') as vcf:
    for _ in range(headerlines):
        next(vcf)
    if filt == True:
        for line in vcf.readlines():
            if line.rsplit()[6] == 'PASS':
                lineout = parser(line, nsamples, maf)
                if lineout != '':
                    for i in range(nsamples):
                        output[i] += lineout[i]
    else:
        for line in vcf.readlines():
            lineout = parser(line, nsamples, maf)
            if lineout != '':
                for i in range(nsamples):
                    output[i] += lineout[i]

# Final printing with PHYLIP standart specifications
print('   {}   {}'.format(nsamples, len(output[0])))
for nsample in range(nsamples):
    lenspaces = 10 - len(samples[nsample])
    print('{}{}{}'.format(samples[nsample], " "*lenspaces, output[nsample]))
