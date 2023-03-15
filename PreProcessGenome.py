import os, sys
import csv as csv
import numpy as np
from Bio.Seq import Seq
import optparse
import logging

codonDict = {'GCT': 0b00000000000000000000000001,
			 'GCC': 0b00000000000000000000000001,
			 'GCA': 0b00000000000000000000000001,
			 'GCG': 0b00000000000000000000000001,
			 'CGT': 0b00000000000000000000000010,
			 'CGC': 0b00000000000000000000000010,
			 'CGA': 0b00000000000000000000000010,
			 'CGG': 0b00000000000000000000000010,
			 'AGA': 0b00000000000000000000000010,
			 'AGG': 0b00000000000000000000000010,
			 'AAT': 0b00000000000000000000000100,
			 'AAC': 0b00000000000000000000000100,
			 'GAT': 0b00000000000000000000001000,
			 'GAC': 0b00000000000000000000001000,
			 'TGT': 0b00000000000000000000010000,
			 'TGC': 0b00000000000000000000010000,
			 'CAA': 0b00000000000000000000100000,
			 'CAG': 0b00000000000000000000100000,
			 'GAA': 0b00000000000000000001000000,
			 'GAG': 0b00000000000000000001000000,
			 'GGT': 0b00000000000000000010000000,
			 'GGC': 0b00000000000000000010000000,
			 'GGA': 0b00000000000000000010000000,
			 'GGG': 0b00000000000000000010000000,
			 'CAT': 0b00000000000000000100000000,
			 'CAC': 0b00000000000000000100000000,
			 'ATT': 0b00000000000000001000000000,
			 'ATC': 0b00000000000000001000000000,
			 'ATA': 0b00000000000000001000000000,
			 'TTA': 0b00000000000000010000000000,
			 'TTG': 0b00000000000000010000000000,
			 'CTT': 0b00000000000000010000000000,
			 'CTC': 0b00000000000000010000000000,
			 'CTA': 0b00000000000000010000000000,
			 'CTG': 0b00000000000000010000000000,
			 'AAA': 0b00000000000000100000000000,
			 'AAG': 0b00000000000000100000000000,
			 'ATG': 0b00000000000001000000000000,
			 'TTT': 0b00000000000010000000000000,
			 'TTC': 0b00000000000010000000000000,
			 'CCT': 0b00000000000100000000000000,
			 'CCC': 0b00000000000100000000000000,
			 'CCA': 0b00000000000100000000000000,
			 'CCG': 0b00000000000100000000000000,
			 'TCT': 0b00000000001000000000000000,
			 'TCC': 0b00000000001000000000000000,
			 'TCA': 0b00000000001000000000000000,
			 'TCG': 0b00000000001000000000000000,
			 'AGT': 0b00000000001000000000000000,
			 'AGC': 0b00000000001000000000000000,
			 'ACT': 0b00000000010000000000000000,
			 'ACC': 0b00000000010000000000000000,
			 'ACA': 0b00000000010000000000000000,
			 'ACG': 0b00000000010000000000000000,
			 'TGG': 0b00000000100000000000000000,
			 'TAT': 0b00000001000000000000000000,
			 'TAC': 0b00000001000000000000000000,
			 'GTT': 0b00000010000000000000000000,
			 'GTC': 0b00000010000000000000000000,
			 'GTA': 0b00000010000000000000000000,
			 'GTG': 0b00000010000000000000000000,
			 'TAA': 0b00000100000000000000000000,
			 'TAG': 0b00000100000000000000000000,
			 'TGA': 0b00000100000000000000000000,
			 'A': 0b00001000000000000000000000,
			 'T': 0b00010000000000000000000000,
			 'C': 0b00100000000000000000000000,
			 'G': 0b01000000000000000000000000,
			 'ERROR': 0b10000000000000000000000000,
			 }

codonAminoDict = {'GCT': 'A',
			 'GCC': 'A',
			 'GCA': 'A',
			 'GCG': 'A',
			 'CGT': 'R',
			 'CGC': 'R',
			 'CGA': 'R',
			 'CGG': 'R',
			 'AGA': 'R',
			 'AGG': 'R',
			 'AAT': 'N',
			 'AAC': 'N',
			 'GAT': 'D',
			 'GAC': 'D',
			 'TGT': 'C',
			 'TGC': 'C',
			 'CAA': 'Q',
			 'CAG': 'Q',
			 'GAA': 'E',
			 'GAG': 'E',
			 'GGT': 'G',
			 'GGC': 'G',
			 'GGA': 'G',
			 'GGG': 'G',
			 'CAT': 'H',
			 'CAC': 'H',
			 'ATT': 'I',
			 'ATC': 'I',
			 'ATA': 'I',
			 'TTA': 'L',
			 'TTG': 'L',
			 'CTT': 'L',
			 'CTC': 'L',
			 'CTA': 'L',
			 'CTG': 'L',
			 'AAA': 'K',
			 'AAG': 'K',
			 'ATG': 'M',
			 'TTT': 'F',
			 'TTC': 'F',
			 'CCT': 'P',
			 'CCC': 'P',
			 'CCA': 'P',
			 'CCG': 'P',
			 'TCT': 'S',
			 'TCC': 'S',
			 'TCA': 'S',
			 'TCG': 'S',
			 'AGT': 'S',
			 'AGC': 'S',
			 'ACT': 'T',
			 'ACC': 'T',
			 'ACA': 'T',
			 'ACG': 'T',
			 'TGG': 'W',
			 'TAT': 'Y',
			 'TAC': 'Y',
			 'GTT': 'V',
			 'GTC': 'V',
			 'GTA': 'V',
			 'GTG': 'V',
			 'TAA': 'Z',
			 'TAG': 'Z',
			 'TGA': 'Z',
			  'ERROR':'X'
			 }

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-i", "--fileFolder", action = "store", type = "string", dest = "fileFolder",
									help = "fileFolder")
parser.add_option("-o", "--outFolder", action = "store", type = "string", dest = "outFolder",
									help = "outFolder")

(options, args) = parser.parse_args()
if options.fileFolder is None:
	sys.stderr.write(prog_base + ": ERROR: missing required command-line argument")
	parser.print_help()
	sys.exit(0)

fileFolder = options.fileFolder
outFolder = options.outFolder

class CDSEntity:
	def __init__(self, region, translation):
		self.__region = region
		self.__region_temp = region
		self.__translation = translation
		if self.__region.startswith('complement'):
			self.__complement = True
			self.__region_temp = self.__region.strip('complement()')
		else:
			self.__complement = False
		self.__comment = ''
		self.__unmatch_percent = 0
		self.__X_percent = 0
		self.__Z_percent = 0
		self.__trim_endZ = False
		self.__refined_region = False
		self.__region_list = []
		if self.__region_temp.startswith('join'):
			region_strings = self.__region_temp.strip('join()').split(',')
			for region_string in region_strings:
				if region_string.startswith('complement'):
					regions = region_string.strip('complement()').split('..')
					if len(regions) == 2:
						self.__region_list.append([regions[0].strip('<>'), regions[1].strip('<>'), True])
					else:
						self.__region_list.append([regions[0].strip('<>'), regions[0].strip('<>'), True])
				else:
					regions = region_string.split('..')
					if len(regions) == 2:
						self.__region_list.append([regions[0].strip('<>'), regions[1].strip('<>')])
					else:
						self.__region_list.append([regions[0].strip('<>'), regions[0].strip('<>')])
		else:
			regions = self.__region_temp.split('..')
			self.__region_list.append([regions[0].strip('<>'), regions[1].strip('<>')])

	@property
	def region(self):
		return self.__region

	@property
	def region_list(self):
		return self.__region_list

	@region_list.setter
	def region_list(self, region_list):
		self.__region_list = region_list

	@property
	def complement(self):
		return self.__complement

	@property
	def translation(self):
		return self.__translation

	@translation.setter
	def translation(self, translation):
		self.__translation = translation

	@property
	def comment(self):
		return self.__comment
	@comment.setter
	def comment(self, comment):
		self.__comment = comment

	@property
	def unmatch_percent(self):
		return self.__unmatch_percent
	@unmatch_percent.setter
	def unmatch_percent(self, unmatch_percent):
		self.__unmatch_percent = round(unmatch_percent, 2)

	@property
	def X_percent(self):
		return self.__X_percent

	@X_percent.setter
	def X_percent(self, X_percent):
		self.__X_percent = round(X_percent, 2)

	@property
	def Z_percent(self):
		return self.__Z_percent

	@Z_percent.setter
	def Z_percent(self, Z_percent):
		self.__Z_percent = round(Z_percent, 2)

	@property
	def trim_endZ(self):
		return self.__trim_endZ
	@trim_endZ.setter
	def trim_endZ(self, trim_endX):
		self.__trim_endZ = trim_endX

	@property
	def refined_region(self):
		return self.__refined_region

	@refined_region.setter
	def refined_region(self, refined_region):
		self.__refined_region = refined_region

class GeneAnnotation:
	def __init__(self, name, length, code_type, shape, date):
		self.__name = name
		self.__length = length
		self.__code_type = code_type
		self.__shape = shape
		self.__date = date
		self.__accession = ''
		self.__version = ''
		self.__CDS = []
		self.__sequence = ''
		self.__author = ''

	@property
	def name(self):
		return self.__name

	@property
	def length(self):
		return self.__length

	@property
	def code_type(self):
		return self.__code_type

	@property
	def shape(self):
		return self.__shape

	@property
	def date(self):
		return self.__date

	@property
	def accession(self):
		return self.__accession

	@accession.setter
	def accession(self, accession):
		self.__accession = accession

	@property
	def version(self):
		return self.__version

	@version.setter
	def version(self, version):
		self.__version = version

	@property
	def authors(self):
		return self.__author

	@authors.setter
	def authors(self, author):
		self.__author = author

	@property
	def CDS(self):
		return self.__CDS

	@CDS.setter
	def CDS(self, cds_entity):
		self.__CDS.append(cds_entity)

	@property
	def sequence(self):
		return self.__sequence

	@sequence.setter
	def sequence(self, sequence):
		self.__sequence = sequence


class Logger:
	def __init__(self, name):
		logging.basicConfig(
			level = logging.DEBUG,
			format = '%(asctime)s  %(filename)s : %(levelname)s  %(message)s',
			datefmt = '%Y-%m-%d %A %H:%M:%S',
			filename = name,
			filemode = 'a'
		)
		console = logging.StreamHandler()
		console.setLevel(logging.INFO)
		formatter = logging.Formatter('%(asctime)s  %(filename)s : %(levelname)s  %(message)s')
		console.setFormatter(formatter)
		# Create an instance
		logging.getLogger().addHandler(console)

	def log_info(self, message):
		logging.info(message)

	def log_warning(self, message):
		logging.warning(message)

	def log_error(self, message):
		logging.error(message)


def encodeSeqSegment(seq, start, end, is_cds):
	seqCode = list()
	if is_cds:
		pos = start - 1
		while pos + 3 <= end:
			seqCode.append(codonDict[seq[pos:pos + 3]])
			pos += 3
	else:
		for pos in range(start - 1, end - 1):
			seqCode.append(codonDict[seq[pos]])
	return seqCode

def encodeSequence(sequence, cdses, length):
	start = 1
	enSeq = []
	for cds in cdses:
		startPos = cds[0]
		endPos = cds[1]
		if startPos > start:
			enSeq = np.append(enSeq, encodeSeqSegment(sequence, start, startPos - 1, False))
		enSeq = np.append(enSeq, encodeSeqSegment(sequence, startPos, endPos, True))
		start = endPos + 1
	if start < length:
		enSeq = np.append(enSeq, encodeSeqSegment(sequence, start, length, False))
	return enSeq


def match_sequence(seqA, seqB):
	equal_start = 0
	equal_end = -1
	replace_start = 0
	replace_end = -1
	opcodes = []
	is_last_equal = False
	replace_count = 0
	equal_count = 0
	max_equal = 0
	min_len = len(seqA)
	if min_len > len(seqB):
		min_len = len(seqB)
	except_chars = ['X', 'Z']

	for i in range(0, min_len):
		if seqA[i] == seqB[i] or seqA[i] in except_chars or seqB[i] in except_chars:
			if not is_last_equal:
				equal_start = i
				if equal_start == replace_end:
					opcodes.append(('replace', replace_start, replace_end))
					replace_count += replace_end - replace_start
			equal_end = i + 1
			is_last_equal = True
		elif seqA[i] != seqB[i]:
			if is_last_equal:
				replace_start = i
				if replace_start == equal_end:
					opcodes.append(('equal', equal_start, equal_end))
					if replace_count > 0:
						equal_count += equal_end - equal_start

					if equal_end - equal_start > max_equal:
						max_equal = equal_end - equal_start
			replace_end = i + 1
			is_last_equal = False

		if i == len(seqA) - 1:
			if is_last_equal:
				opcodes.append(('equal', equal_start, equal_end))
				equal_count += equal_end - equal_start

				if equal_end - equal_start > max_equal:
					max_equal = equal_end - equal_start
			else:
				opcodes.append(('replace', replace_start, replace_end))
				replace_count += replace_end - replace_start

	if min_len < len(seqA):
		opcodes.append(('delete', min_len, len(seqA)))
	elif min_len < len(seqB):
		opcodes.append(('insert', min_len, len(seqB)))

	replace_exceed = replace_count > equal_count and equal_count > 0
	if not replace_exceed: # ajust replace_exceed
		if max_equal < len(seqA) / 2:
			replace_exceed = True
	return opcodes, replace_exceed


def refine_cds(cds_entity, codon_start, i):
	length = 0
	if cds_entity.complement:
		for region_index in range(len(cds_entity.region_list) - 1, -1, -1):
			if length + int(cds_entity.region_list[region_index][1]) - int(
					cds_entity.region_list[region_index][0]) + 1 > codon_start:
				region_start = cds_entity.region_list[region_index][0]
				cds_entity.region_list[region_index][0] = int(cds_entity.region_list[region_index][1]) - codon_start + i - length + 1
				cds_entity.region_list.insert(region_index,
											  [region_start, int(cds_entity.region_list[region_index][1]) - codon_start - length])
				if int(cds_entity.region_list[region_index + 1][1]) < cds_entity.region_list[region_index + 1][0]:
					cds_entity.region_list.pop(region_index + 1)
				break
			else:
				length += int(cds_entity.region_list[region_index][1]) - int(
					cds_entity.region_list[region_index][0]) + 1
	else:
		for region_index in range(0, len(cds_entity.region_list)):
			if length + int(cds_entity.region_list[region_index][1]) - int(cds_entity.region_list[region_index][0]) + 1 > codon_start:
				region_end = cds_entity.region_list[region_index][1]
				cds_entity.region_list[region_index][1] = codon_start - i - length + int(cds_entity.region_list[region_index][0]) - 1
				cds_entity.region_list.insert(region_index + 1, [codon_start - length + int(cds_entity.region_list[region_index][0]), region_end])
				if cds_entity.region_list[region_index][1] < int(cds_entity.region_list[region_index][0]):
					cds_entity.region_list.pop(region_index)
				break
			else:
				length += int(cds_entity.region_list[region_index][1]) - int(cds_entity.region_list[region_index][0]) + 1
	cds_entity.refined_region = True

def match_without_XZ(translation, trans_string):
	match = True
	min_len = min(len(translation), len(trans_string))
	for index in range(0, min_len):
		if translation[index] != 'X' and trans_string[index] != 'X'\
			and trans_string[index] != 'Z' and translation[index] != trans_string[index]:
			match = False
			break
	return match, len(translation) - len(trans_string)

def adjust_start(gene_annotation, cds_entity, translation, cds_seq, start):
	for i in range(-5, 6):
		if i == 0:
			continue
		codon_start = start + i
		if i < 0 and not cds_entity.complement and codon_start + int(cds_entity.region_list[0][0]) - 1 < 0:
			continue
		elif i < 0 and cds_entity.complement and int(cds_entity.region_list[-1][1]) - codon_start > len(
				gene_annotation.sequence):
			continue
		elif i < 0 and start == 0:
			continue

		precheck_len = 4
		check_string = ''
		for index in range(0, precheck_len):
			amino_start = codon_start + (index * 3)
			if amino_start + 3 > len(cds_seq):
				precheck_len = index
				break
			if codonAminoDict.__contains__(cds_seq[amino_start:amino_start + 3]):
				check_string += codonAminoDict[cds_seq[amino_start:amino_start + 3]]
			else:
				check_string += 'X'

		match, diff = match_without_XZ(translation[:precheck_len], check_string)
		if match:
			return i
	return 0

def handle_region_longer(cds_entity, gap):
	if cds_entity.complement:
		cds_entity.region_list[0][0] = int(cds_entity.region_list[0][0]) + gap
		cds_entity.comment += 'cds region has error, trim {0} at the start;'.format(gap)
	else:
		cds_entity.region_list[-1][1] = int(cds_entity.region_list[-1][1]) + gap
		cds_entity.comment += 'cds region has error, trim {0} at the end;'.format(gap)
	cds_entity.refined_region = True

def handle_translation_longer(cds_entity, gap, encodeTrans, start, end):
	if cds_entity.complement:
		cached_start = int(cds_entity.region_list[0][0])
		if cached_start - gap <= 0:  # if sequence is not complete for translation, trim the translation and sequence
			cds_entity.translation = cds_entity.translation[0:-1]
			cds_entity.region_list[0][0] = cached_start + 3 - gap
			cds_entity.comment += 'sequence is not long enough, trim the first in cds [{0}..{1}] translation, and change to {2};'.format(
				cds_entity.region_list[0][0], cached_start, cds_entity.region_list[0][0])
		else:
			cds_entity.region_list[0][0] = cached_start - 3 + gap
			cds_entity.comment += 'cds region has error, extend {0} at the start;'.format(gap)
			encodeTrans.extend(translation[start: end])
	else:
		cached_end = int(cds_entity.region_list[-1][1])
		if cached_end + gap > len(gene_annotation.sequence):
			cds_entity.translation = cds_entity.translation[0:-1]
			cds_entity.comment += 'sequence is not long enough, trim the cds [{0}..{1}] translation, and change to {2};'.format(
				cds_entity.region_list[-1][0], cached_end, cds_entity.region_list[-1][1])
		else:
			cds_entity.region_list[-1][1] = cached_end + gap
			cds_entity.comment += 'cds region has error, extend {0} at the end;'.format(gap)
			logger.log_info('{0} cds {1} region has error, extend {2} at the end'.format(
				gene_annotation.name, cds_entity.region_list[-1], gap))
			encodeTrans.extend(translation[start: end])
	cds_entity.refined_region = True

def refine_translation(translation, trans_string, cds_seq, gene_annotation, cds_entity, start = 0, adjust = 0):
	encodeTrans = []
	opcodes, replace_exceed = match_sequence(translation, trans_string)

	refined = False
	for opcode_index in range(0, len(opcodes)):
		opcode = opcodes[opcode_index]
		if opcode[0] == 'equal':
			encodeTrans.extend(translation[opcode[1]: opcode[2]])
		elif opcode[0] == 'replace':
			# if not replace length exceeds equal length, or the last replace, or the following equal longer than replace, or the replace length=1, and not the first is replace
			if len(opcodes) > 1 and not replace_exceed and (opcode_index == len(opcodes) - 1\
					or opcode[2] - opcode[1] < opcodes[opcode_index + 1][2] - opcodes[opcode_index + 1][1])\
					or opcode[2] - opcode[1] == 1 and opcode[2] - opcode[1] < opcodes[opcode_index + 1][2] - opcodes[opcode_index + 1][1]:
				encodeTrans.extend(translation[opcode[1]: opcode[2]])
			else:
				adjust = adjust_start(gene_annotation, cds_entity, translation[opcode[1]:], cds_seq,
												  opcode[1] * 3 + start)
				if adjust == 0:
					encodeTrans.extend(translation[opcode[1]: opcode[2]])
					continue
				codon_start = opcode[1] * 3 + adjust + start

				refine_cds(cds_entity, codon_start, adjust)

				refined_trans = []
				for codon_index in range(codon_start, len(cds_seq) - 2, 3):
					if codonAminoDict.__contains__(cds_seq[codon_index:codon_index + 3]):
						refined_trans.extend(codonAminoDict[cds_seq[codon_index:codon_index + 3]])
					else:
						refined_trans.extend('X')  # refine to ERROR codon

				new_translation = translation[opcode[1]:]
				new_trans_string = ''.join(refined_trans).rstrip('Z')
				match, diff = match_without_XZ(new_translation, new_trans_string)
				if match:
					encodeTrans.extend(refined_trans)
					if diff > 0:
						gap = start + len(translation) * 3 - len(cds_seq) + adjust
						handle_translation_longer(cds_entity, gap, encodeTrans, opcode[1], opcode[2])
					elif diff < 0:
						gap = start + len(translation) * 3 - len(cds_seq) + adjust
						handle_region_longer(cds_entity, gap)
					refined = True
				else:
					logger.log_info('recursive refine_translation {0}, index {1}'.format(opcodes, opcode_index))
					refine_translation(new_translation, new_trans_string, cds_seq, gene_annotation,
									   cds_entity, codon_start, adjust)
					refined = True
				break
		elif opcode[0] == 'insert':
			gap = start + len(translation) * 3 - len(cds_seq) + adjust
			handle_region_longer(cds_entity, gap)
		elif opcode[0] == 'delete':
			if start + len(translation) * 3 > len(cds_seq):
				gap = start + len(translation) * 3 - len(cds_seq) + adjust
				handle_translation_longer(cds_entity, gap, encodeTrans, opcode[1], opcode[2])

	return encodeTrans

def translate_sequence(gene_annotation, cds_entity, record_comment):
	encodeTrans = list()
	cds_seq = ''
	unmatch_count = 0
	if cds_entity.complement:
		for cds in reversed(cds_entity.region_list):
			start = int(cds[0])
			end = int(cds[1])
			cds_seq += Seq(gene_annotation.sequence[start - 1:end]).reverse_complement()
	else:
		for cds in cds_entity.region_list:
			start = int(cds[0])
			end = int(cds[1])
			if len(cds) > 2:
				# complement sub region like join(153006..153821,complement(82031..82618),complement(80481..81305))
				cds_seq += Seq(gene_annotation.sequence[start - 1:end]).reverse_complement()
			else:
				cds_seq += gene_annotation.sequence[start - 1:end]
	pos = 0
	while pos + 3 <= len(cds_seq):
		condon = cds_seq[pos:pos + 3]
		pos_trans = int(pos / 3)
		if codonAminoDict.__contains__(condon):
			encoded_protein = codonAminoDict[condon]
			if record_comment:
				if len(cds_entity.translation) > pos_trans and encoded_protein != cds_entity.translation[pos_trans] \
						or len(cds_entity.translation) <= pos_trans:
					unmatch_count += 1
			encodeTrans.append(encoded_protein)
		else:
			if pos_trans >= len(cds_entity.translation):
				encodeTrans.append('Z')  # only few sequence (2 in 23W) has such error, so refine to end codon
				if record_comment:
					cds_entity.comment += '{0} is translated to Z at the end;'.format(condon)
					unmatch_count += 1
			else:
				# uses translation code
				encodeTrans.append(cds_entity.translation[pos_trans])
				if record_comment:
					cds_entity.comment += '{0} is translated to {1};'.format(condon, cds_entity.translation[pos_trans])
					unmatch_count += 1
		pos += 3

	if record_comment and len(cds_entity.translation) > int(pos / 3):
		unmatch_count += len(cds_entity.translation) - int(pos / 3)
	if record_comment and unmatch_count > 0:
		cds_entity.unmatch_percent = unmatch_count / len(cds_entity.translation) * 100

	return encodeTrans, cds_seq


def checkTrans(gene_annotation):
	for cds_entity in gene_annotation.CDS:
		encodeTrans, cds_seq = translate_sequence(gene_annotation, cds_entity, False)

		trans_string = ''
		if encodeTrans[-1] == 'Z':
			trans_string = ''.join(encodeTrans[0:-1])
		else:
			trans_string = ''.join(encodeTrans)
		if trans_string != cds_entity.translation:
			if abs(len(trans_string) - len(cds_entity.translation)) > 5:
				logger.log_warning(
					'{0} length {1} CDS complement:{2} [{3}] difference of length > 5 translation len: {4} encode len: {5}, will not do refine'
					.format(gene_annotation.name, gene_annotation.length, cds_entity.complement,
							cds_entity.region_list, len(cds_entity.translation), len(trans_string),
							))
				continue
			refine_translation(cds_entity.translation, trans_string, cds_seq, gene_annotation, cds_entity)

		encodeTrans, seq2 = translate_sequence(gene_annotation, cds_entity, True)
		if encodeTrans[-1] == 'Z':
			trans_string = ''.join(encodeTrans[0:-1])
			cds_entity.trim_endZ = True
		else:
			trans_string = ''.join(encodeTrans)

		if trans_string != cds_entity.translation:
			if trans_string.find('Z') > 0:
				Z_percent = np.char.count(trans_string, 'Z') / len(trans_string) * 100
				cds_entity.Z_percent = Z_percent
				if cds_entity.Z_percent > 5:
					logger.log_info('{0} length {1} CDS complement:{2} [{3}] Z percent {4}%'.format(gene_annotation.name, gene_annotation.length, cds_entity.complement,
							cds_entity.region_list, cds_entity.Z_percent))
			if (cds_entity.unmatch_percent > 10):
				logger.log_info(cds_entity.translation)
				logger.log_info(trans_string)
				logger.log_info(
					'{0} length {1} CDS complement:{2} [{3}] unmatch percent={4}%, translation len: {5} encode len: {6}'
					.format(gene_annotation.name, gene_annotation.length, cds_entity.complement,
							cds_entity.region_list, cds_entity.unmatch_percent, len(cds_entity.translation), len(trans_string),
							))
		X_percent = np.char.count(cds_entity.translation, 'X') / len(cds_entity.translation) * 100
		cds_entity.X_percent = X_percent
		if X_percent > 5:
			logger.log_info(
				'{0} length {1} CDS complement:{2} [{3}] X percent={4}%, translation len: {5} encode len: {6}'
				.format(gene_annotation.name, gene_annotation.length, cds_entity.complement,
						cds_entity.region_list, cds_entity.X_percent, len(cds_entity.translation), len(trans_string),
						))

gene_sequences_dic = dict()
file_cached = dict()
def read_sequence(gbff_filename, id, contig_region):
	if id not in gene_sequences_dic:
		filename_sections = gbff_filename.split('.')
		if len(filename_sections) == 4:
			fna_filename = '.'.join(filename_sections[:3]) + '.fna'
			if fna_filename not in file_cached:
				if len(file_cached) > 0 and gbff_filename not in file_cached.values():
					file_cached.clear()
					gene_sequences_dic.clear()
					logger.log_info('clean file and sequence cache')

				hit = False
				i = 1
				fna_filepath = os.path.join(fileFolder, fna_filename)
				while not hit:
					if os.path.exists(fna_filepath):
						logger.log_info('start to read sequences from file: ' + fna_filename)
						with open(fna_filepath, 'r') as fna_file:
							line = fna_file.readline()
							while line != '':
								if line.startswith('>'):
									sequence = ''
									infos = line.lstrip('>').split()
									line = fna_file.readline()
									while not line.startswith('>') and line != '':
										sequence += line.strip()
										line = fna_file.readline()
									if infos[0] == id:
										hit = True
									gene_sequences_dic[infos[0]] = sequence
						logger.log_info('finished to read file: ' + fna_filename)
					if not hit:
						fna_filename = '{0}.{1}.{2}.fna'.format('.'.join(filename_sections[:2]), i, filename_sections[2])
						fna_filepath = os.path.join(fileFolder, fna_filename)
						i += 1
				file_cached[fna_filename] = gbff_filename

	regions = contig_region.rstrip(')').split(':')[1].split('..')
	return gene_sequences_dic[id][int(regions[0]) - 1: int(regions[1])]


def save_gene_cds(gene_annotation, no_cds_gene_writer):
	if len(gene_annotation.CDS) == 0:
		no_cds_gene_writer.writerow([gene_annotation.name, gene_annotation.length, gene_annotation.code_type,
						  gene_annotation.shape, gene_annotation.date, gene_annotation.accession,
						  gene_annotation.version, gene_annotation.authors, gene_annotation.sequence])
		logger.log_warning('{0} has no quality CDS'.format(gene_annotation.name))
		return
	checkTrans(gene_annotation)
	gene_writer.writerow([gene_annotation.name, gene_annotation.length, gene_annotation.code_type,
						  gene_annotation.shape, gene_annotation.date, gene_annotation.accession,
						  gene_annotation.version, gene_annotation.authors, gene_annotation.sequence])
	for cds_entity in gene_annotation.CDS:
		cds_writer.writerow([gene_annotation.accession, cds_entity.region, cds_entity.region_list,
							 cds_entity.refined_region,
							 cds_entity.complement, cds_entity.unmatch_percent, cds_entity.X_percent,
							 cds_entity.Z_percent, cds_entity.trim_endZ, cds_entity.comment,
							 cds_entity.translation])


logger = Logger('process.log')
files = os.listdir(fileFolder)
for file_name in files:
	fname, ext = os.path.splitext(file_name)
	if (ext != '.gbff'):
		continue
	logger.log_info('start to process file ' + file_name)
	file_path = os.path.join(fileFolder, file_name)
	gene_file_path = os.path.join(outFolder, fname + '_gene_file.csv')
	cds_file_path = os.path.join(outFolder, fname + '_cds_file.csv')
	no_cds_gene_file_path = os.path.join(outFolder, fname + '_no_cds_gene.csv')
	with open(gene_file_path, 'w', encoding='utf-8', newline='') as gene_file,\
			open(cds_file_path, 'w', encoding='utf-8', newline='') as cds_file, \
			open(no_cds_gene_file_path, 'w', encoding='utf-8', newline='') as no_cds_gene_file, \
			open(file_path, 'r') as gbff_file:
		gene_writer = csv.writer(gene_file)
		gene_writer.writerow(['name', 'length', 'type', 'shape', 'submit date', 'accession', 'version', 'authors' ,'sequence'])
		no_cds_gene_writer = csv.writer(no_cds_gene_file)
		no_cds_gene_writer.writerow(
			['name', 'length', 'type', 'shape', 'submit date', 'accession', 'version', 'authors','sequence'])
		cds_writer = csv.writer(cds_file)
		cds_writer.writerow(
			['accession', 'cds region', 'region list', 'refined region', 'complement', 'unmatch percent', 'X percent in trans',
			 'Z percent in encode', 'end with terminate codon', 'comment', 'translation'])

		gene_annotation = None
		line = gbff_file.readline()
		while line != '':
			if line.startswith('LOCUS'):
				locus_strings = line.split()
				gene_annotation = GeneAnnotation(locus_strings[1], locus_strings[2], locus_strings[4], locus_strings[5], locus_strings[7])
			elif line.startswith('ACCESSION'):
				gene_annotation.accession = line.split()[1]
			elif line.startswith('VERSION'):
				gene_annotation.version = line.split()[1]
			elif line.startswith('  AUTHORS'):
				gene_annotation.authors = line.split()[1]
			elif line.startswith('     CDS'):
				region = line.split()[1]
				next_line = gbff_file.readline().strip()
				while not next_line.startswith('/'):
					region += next_line
					next_line = gbff_file.readline().strip()
				has_translation = True
				while not next_line.startswith('/translation='):
					next_line = gbff_file.readline().strip()
					if next_line.startswith('CONTIG') or next_line.startswith('gene            '):
						has_translation = False
						logger.log_warning('{0} CDS {1} has no translation'.format(gene_annotation.name, region))
						break
				if has_translation:
					translation = next_line.lstrip('/translation="')
					while not next_line.endswith('"'):
						next_line = gbff_file.readline().strip()
						translation += next_line
					translation = translation.strip('"')
					gene_annotation.CDS = CDSEntity(region, translation)
				elif next_line.startswith('CONTIG'):
					gene_annotation.sequence = read_sequence(file_name, gene_annotation.version, next_line.split()[1])
					save_gene_cds(gene_annotation, no_cds_gene_writer)
			elif line.startswith('CONTIG'):
				gene_annotation.sequence = read_sequence(file_name, gene_annotation.version, line.split()[1])
				save_gene_cds(gene_annotation, no_cds_gene_writer)
			elif line.startswith('ORIGIN'):
				sequence = ''
				while True:
					line = gbff_file.readline().strip()
					if line == '//':
						break
					segments = line.split()
					segments.pop(0)
					sequence += ''.join(segments)
				gene_annotation.sequence = sequence.upper()
				save_gene_cds(gene_annotation, no_cds_gene_writer)

			line = gbff_file.readline()
	logger.log_info('finished to process file ' + file_name)



