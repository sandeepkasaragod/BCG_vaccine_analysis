from argparse import ArgumentParser
from os.path import join as join
import os

def convert(args):
	root_seq = read_snp_caller(args.snp_caller_file, args.genome_name, args.genome_version)
	dicts = {}
	
	os.makedirs(args.output,exist_ok=True)

	write_file = open(args.output + '/' + args.genome_name + '_placed.infile.fasta', 'w')
	write_file.write('>' + root_seq[0] + '\n' + root_seq[1].rstrip() + '\n')

	for i in open(args.input + '_placed.infile'):
		split_i = i.split(' ')
		if split_i[0] not in dicts:
			dicts[split_i[0]] = 1
			write_file.write('>' + split_i[0] + '\n' + split_i[2].rstrip() + '\n')
	write_file.close()

def read_snp_caller(infile, genome_name, version):
	for i in open(infile):
		split_i = i.split(' ')
		if genome_name + '.' + version in split_i[0] in split_i[0]:
			return split_i[0], split_i[2]
			break
	
if __name__ =="__main__":
	parser = ArgumentParser(description='convert to fasta')
	parser.add_argument('-i', '--input', help='input file', required=True)
	parser.add_argument('-o', '--output', help='output fasta file', required=True)
	parser.add_argument('-s', '--snp_caller_file', help='output from snp caller', required=True)
	parser.add_argument('-g', '--genome_name', help='genome name', required=True)
	parser.add_argument('-v', '--genome_version', help='genome version', required=True)
	args = parser.parse_args()
	convert(args)
