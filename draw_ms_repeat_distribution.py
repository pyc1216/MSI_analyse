#!/usr/bin/python
#coding=utf-8:
#########################################
#VERSION 1.2
#AUTHOR yu.pu
#DATA 03-18-2020
#读取 *msi.pcr.result_dis, 画出5个pcr位点的repeat times (归一化后)的分布
'''
卡方检验要求：最好是大样本数据。一般每个个案最好出现一次，四分之一的个案至少出现五次。如果数据不符合要求，就要应用校正卡方。
总结
1.所有的理论数T≥5并且总样本量n≥40,用Pearson卡方进行检验.
2.如果理论数T＜5但T≥1,并且n≥40,用连续性校正的卡方进行检验.
3.如果有理论数T＜1或n＜40,则用Fisher’s检验.

作者：thinkando
链接：https://www.jianshu.com/p/f0e1b0100e59
来源：简书

'''
#########################################

from __future__ import print_function
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency
import math


## key=pcr msi site name, value=position in genome(hg19)
ALL_MS = {
'MONO27': 39536689, 
'BAT26_5': 47641559, 
'NR24': 95849361, 
'BAT25': 55598211, 
'NR21': 23652346
}

def get_pass_index(dis):
	#2.remove outliers (mean ± 3 * std)
	values = []
	index = []
	for i, cnt in enumerate(dis):
		k = i + 1
		values += [k] * cnt
	mean = np.mean(values)
	std = np.std(values)
	sds = 3
	min = math.floor(mean - sds * std)
	max = math.ceil(mean + sds * std)
	for i, cnt in enumerate(dis):
		k = i + 1
		if k <= min or k >= max:
			continue
		index.append(i)
	return index



def cal_chi2_pvalue(normal_dis, tumor_dis):
	'''
	input distribution:
	  normal_dis = [0, 0, 4, 10, 11, 8, 25, 36, 64, 96, 94, 91, 39, 14, 9, 1, 1, 0, 0]
	  tumor_dis =  [0, 0, 0, 10, 17, 35, 58, 105, 129, 191, 253, 216, 87, 24, 9, 4, 4, 4, 0]
	filter:
	  1.remove 0 in head/tail
	  2.remove outliers (mean ± 3 * std)
	  3.normalize(expected total reads = 500)
	  4.remove chi2 test's expected value < 5
	output distribution:
	  normal_dis_filter2 = [7.62, 10.97, 19.5, 29.25, 28.64, 27.73, 11.88]
	  tumor_dis_filter2 = [17.67, 31.99, 39.31, 58.2, 77.09, 65.81, 26.51]
	'''
	#1.
	index_both_zero = []
	for i, count in enumerate(zip(normal_dis, tumor_dis)):
		if count == (0,0):
			index_both_zero.append(i)
	print('[INFO] normal_dis = ', end='', file=sys.stderr)
	print([x for i, x in enumerate(normal_dis) if i not in index_both_zero], file=sys.stderr)
	print('[INFO] tumor_dis = ', end='', file=sys.stderr)
	print([x for i, x in enumerate(tumor_dis) if i not in index_both_zero], file=sys.stderr)
	
	#2.
	pass_index_n = get_pass_index(normal_dis)
	pass_index_t = get_pass_index(tumor_dis)
	pass_index = list((set(pass_index_n) | set(pass_index_t)) - set(index_both_zero))
	normal_dis_filter = [x for i, x in enumerate(normal_dis) if i in pass_index]
	tumor_dis_filter = [x for i, x in enumerate(tumor_dis) if i in pass_index]
	#3.
	n_reads = sum(normal_dis_filter)
	t_reads = sum(tumor_dis_filter)
	need_reads = 500.0 #magic number
	normal_dis_scale = [x * need_reads / (n_reads + t_reads) for x in normal_dis_filter]
	tumor_dis_scale  = [x * need_reads / (n_reads + t_reads) for x in tumor_dis_filter]
	#4.
	data = np.array([normal_dis_scale, tumor_dis_scale])
	chi2_v, p_value, dof, exp_arr  = chi2_contingency(data)
	normal_exp = exp_arr[0, :]
	tumor_exp = exp_arr[1, :]
	pass_index2 = [all(x) for x in zip(normal_exp >= 5, tumor_exp >= 5)]
	normal_dis_filter2 = data[0,:][pass_index2]
	tumor_dis_filter2 = data[1,:][pass_index2]
	print('[INFO] normal_dis_filter2 = ', end='', file=sys.stderr)
	print([round(x, 2) for x in normal_dis_filter2], file=sys.stderr)
	print('[INFO] tumor_dis_filter2 = ', end='', file=sys.stderr)
	print([round(x, 2) for x in tumor_dis_filter2], file=sys.stderr)
	
	data2 = np.array([normal_dis_filter2, tumor_dis_filter2])

	try:
		chi2_v2, p_value2, dof2, exp_arr2  = chi2_contingency(data2)
	except ValueError:
		chi2_v2, p_value2, dof2, exp_arr2  = chi2_contingency(data)
		print('[WARRING] Filter #4 not applied, p_value({}) not reliable!'.format(p_value2), file=sys.stderr)
	return p_value2


def draw(dict_ms, out_fig):
	'''
	dict_ms = {'MONO27': [[0, 0.5, 0.5, 0.5, 0.2, 1...],[0, 0.3, 0.8, 1.0, 0.5, 0.2, 0...]], ...}
	'''
	#print(dict_ms)
	total_num = len(dict_ms)
	plt.legend(bbox_to_anchor=(1.0, 1), loc=1, borderaxespad=0.)
	#Calculate adjusted p-value
	p_values = []
	d_ranks = {}
	for ms in dict_ms:
		print('[INFO] ms={},pos={}'.format(ms, ALL_MS[ms]), file=sys.stderr)
		normal_dis = dict_ms[ms][0]
		tumor_dis = dict_ms[ms][1]
		p_value = cal_chi2_pvalue(normal_dis, tumor_dis)
		p_values.append((ms, p_value))
		print('*' * 100, file=sys.stderr)
	p_values.sort(key=lambda x: x[1])
	for i, (ms, p_value) in enumerate(p_values):
		rank = i + 1
		q_value = min(p_value * len(dict_ms) / rank, 0.99)
		d_ranks[ms] = (p_value, rank, q_value)
	
	cnt_msi_h = 0
	for ms, (p_value, rank, q_value) in d_ranks.items():
		if q_value <= 0.05:
			cnt_msi_h += 1
	per_msi_h = cnt_msi_h * 1.0 / len(dict_ms)
	
	for i, ms in enumerate(dict_ms):
		normal_dis = dict_ms[ms][0]
		tumor_dis = dict_ms[ms][1]
		p_value, rank, q_value = d_ranks[ms]
		print('p_value=', p_value, '\tq_value=', q_value, file=sys.stderr)
		#*Important*
		scale_func = lambda l:[x * 50.0 / sum(l) for x in l]
		normal_dis_p = scale_func(normal_dis)
		tumor_dis_p = scale_func(tumor_dis)
		fig_rank = i + 1
		axes = plt.subplot(total_num, 1, fig_rank)
		y_max = max(normal_dis_p + tumor_dis_p)
		len_x = max(len(normal_dis_p), len(tumor_dis_p)) - 40 ## usually is equal to 100 - 40 
		plt.plot(range(1, len_x+1), normal_dis_p[:len_x], 'b', label='normal')
		plt.plot(range(1, len_x+1), tumor_dis_p[:len_x], 'r', label='tumor')
		plt.annotate(s=ms, xy=(40, 0.6 * y_max), fontweight='bold')
		plt.annotate(s='p=' + str(q_value), xy=(40,0.3 * y_max))
		if i != (total_num - 1):
			axes.set_xticks([])
		else:
			plt.xlabel('Microsates Repeat Times')
		if i == total_num / 2: #draw y-label in middle of subplot.
			plt.ylabel('Suppot Reads Number(normalization)')
		if i == 0: #draw legend in top-left
			plt.legend(bbox_to_anchor=(1.0, 1), loc=1, borderaxespad=0.)

	print('{:.2%}'.format(per_msi_h), file=sys.stdout)
	if out_fig:
		if out_fig != 'silent':
			plt.savefig(out_fig)
		else:
			return
	plt.show()
	


def load_dis_file(dis_file):
	dict_ms = {}
	
	#msisensor has two dis_file, split by '\t' or ' '.
	with open(dis_file) as f:
		fst_line = f.readline()
		split_by_tab = True if fst_line.count('\t') > fst_line.count(' ') else False 
	
	with open(dis_file) as f:
		
		lines = f.readlines()
		
		i = 0
		while True:
			if i >= len(lines):
				break
			line = lines[i]
		#for i, line in enumerate(lines):
			for ms in ALL_MS:
				if str(ALL_MS[ms]) in line:
					print('[INFO] Loading site(ms={},pos={})...'.format(ms, ALL_MS[ms]), file=sys.stderr)
					if split_by_tab:
						raw_normal_dis = [int(x.strip()) for x in lines[i + 0].split('\t')[7:]]
						raw_tumor_dis  = [int(x.strip()) for x in lines[i + 1].split('\t')[7:]]
						i += 1
					else:
						raw_normal_dis = [int(x) for x in lines[i + 1].lstrip('N:').strip().split()]
						raw_tumor_dis  = [int(x) for x in lines[i + 2].lstrip('T:').strip().split()]
					if sum(raw_normal_dis) == 0:
						normal_dis = [0] * len(raw_normal_dis)
						print('[WARRING] Zero Depth(normal): ms={}'.format(ms), file=sys.stderr)
					else:
						normal_dis = raw_normal_dis
						if sum(raw_normal_dis) < 100:
							print('[WARRING] Low Depth(normal): ms={}, depth={}'.format(ms, str((sum(raw_normal_dis)))), file=sys.stderr)
					if sum(raw_tumor_dis) == 0:
						tumor_dis = [0] * len(raw_tumor_dis)
						print('[WARRING] Zero Depth(tumor): ms={}'.format(ms), file=sys.stderr)
					else:
						tumor_dis = raw_tumor_dis
						if sum(raw_tumor_dis) < 100:
							print('[WARRING] Low Depth(tumor): ms={}, depth={}'.format(ms, str((sum(raw_tumor_dis)))), file=sys.stderr)
					dict_ms[ms] = [normal_dis, tumor_dis]
			i += 1
	return dict_ms


def main():
	dict_ms = load_dis_file(dis_file)
	#print(dict_ms)
	draw(dict_ms, out_fig)


if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) == 1:
		dis_file = args[0]
		out_fig = None
	elif len(args) == 2:
		dis_file, out_fig = args
	else:
		print('usage: {} msi/sample.msi.pcr.result_dis [out_fig]'.format(sys.argv[0]))
		exit(1)
	main()
