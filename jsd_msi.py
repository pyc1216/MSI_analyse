import os
import sys
import numpy as np
from math import log, sqrt
from collections import OrderedDict


def get_entropy(p, q):
	p = [i * 1.0 / sum(p) for i in p]
	q = [i * 1.0 / sum(q) for i in q]
	#if len(p) != 1:
		#p, q = zip(*filter(lambda x: x[0]!=0 or x[1]!=0, zip(p, q)))
	try: #remove both 0 col
		p, q = zip(*filter(lambda x: x[0]!=0 or x[1]!=0, zip(p, q)))
	except ValueError:
		print(p, q)
	p = p + np.spacing(1)
	q = q + np.spacing(1)
	e_p = sum([- i * log(i, 2) for i in p])
	e_q = sum([- j * log(j, 2) for j in q])
	return (e_p, e_q)
	
	
#读取mantis结果

def parse_count(infile):
	d = OrderedDict()
	with open(infile) as f:
		for line in f:
			if line.startswith('Locus'):
				continue
			locus, repeat, normal, tumor = line.strip().split('\t')
			#print(locus, repeat, normal, tumor)
			#break
			repeat, normal, tumor = [int(x) for x in [repeat, normal, tumor]]
			if locus not in d:
				d[locus] = {'normal': [normal], 'tumor': [tumor], 'repeat': [repeat]}
			else:
				d[locus]['normal'].append(normal)
				d[locus]['tumor'].append(tumor)
				d[locus]['repeat'].append(repeat)
	return d


KLD = lambda p, q: sum([_p * log(_p, 2)- _p * log(_q, 2) for (_p, _q) in zip(p, q)])
def JSD(p, q):
	bug = False
	p = [i * 1.0 / sum(p) for i in p]
	q = [i * 1.0 / sum(q) for i in q]
	#if len(p) != 1:
		#p, q = zip(*filter(lambda x: x[0]!=0 or x[1]!=0, zip(p, q)))
	try:
		p, q = zip(*filter(lambda x: x[0]!=0 or x[1]!=0, zip(p, q)))
	except ValueError:
		bug = True
		#print(p, q)
	m = [0.5 * (_p + _q) for _p, _q in zip(p, q)]
	p = p + np.spacing(1)
	q = q + np.spacing(1)
	m = m + np.spacing(1)
	return (0.5 * KLD(p, m) + 0.5 * KLD(q, m), bug)


def cal_mean_jsd(infile, outfile):
	d_count = parse_count(infile)
	sum_jsd = 0
	cnt_filter_by_entropy = 0
	fw = open(outfile, 'w')
	fw.write('Locus\tJSD\tNormal_Entropy\tTumor_Entropy\n')
	for locus in d_count:
		normals = d_count[locus]['normal']
		tumors = d_count[locus]['tumor']
		jsd, bug = JSD(normals, tumors)
		jsd = round(jsd, 4)
		if bug:
			print('[ERROR] {}'.format(locus), file=sys.stderr)
			print(normals, tumors, file=sys.stderr)
		sqrt_jsd = sqrt(jsd)
		e_normal, e_tumor = get_entropy(normals, tumors)
		#print(e_normal, e_tumor)
		fw.write('{}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(locus, jsd, e_normal, e_tumor))
		if e_normal > e_tumor:
			cnt_filter_by_entropy += 1
			continue
		sum_jsd += sqrt_jsd
	fw.close()
	return (round(sum_jsd / len(d_count), 4), cnt_filter_by_entropy, len(d_count))
	

def main(infile, outfile):
	outdir = os.path.dirname(outfile)
	if outdir and not os.path.exists(outdir):
		os.system('mkdir -p {}'.format(outdir))
	outfile_status = outfile + '.status'
	mean_jsd, cnt_filter_by_entropy, cnt_total = cal_mean_jsd(infile, outfile)
	threshold = 0.12
	status = 'Unstable' if mean_jsd >= threshold else 'Stable'
	with open(outfile_status, 'w') as fw:
		fw.write('Average Metric Value (Abbr)\tValue\tThreshold\tStatus\n')
		fw.write('Jensen-Shannon divergence (JSD)\t{}\t{}\t{}\n'.format(mean_jsd, threshold, status))
		
if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) != 2:
		print('Usage: python {} infile outfile'.format(sys.argv[0]), file=sys.stderr)
		exit(1)
	infile, outfile = args
	main(infile, outfile)


