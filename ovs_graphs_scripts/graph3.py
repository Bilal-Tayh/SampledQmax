#!/usr/bin/python
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
from math import log, ceil
import numpy as np
import sys
import itertools
import math
from matplotlib.legend_handler import HandlerErrorbar
from matplotlib import container
import copy



G = sys.argv[1]
sizes=map(lambda s : 10**int(s), sys.argv[2:])
ss = len(sizes)

ocs = cycler(color=['y']*ss) + cycler(marker=['d']*ss)# + cycler(linestyle=['-.']*ss)
hcs = cycler(color=['b']*ss) + cycler(marker=['o']*ss)# + cycler(linestyle=['-']*ss)
scs = cycler(color=['g']*ss) + cycler(marker=['v']*ss)# + cycler(linestyle=['--']*ss)
qcs = cycler(color=['r']*ss) + cycler(marker=['*']*ss)# + cycler(linestyle=[':']*ss)
las_qcs = cycler(color=['b']*ss) + cycler(marker=['p']*ss)# + cycler(linestyle=[':']*ss)
las9_qcs = cycler(color=['g']*ss) + cycler(marker=['P']*ss)# + cycler(linestyle=[':']*ss)

if G == "CPU":
	df_org = pd.read_csv('../results/timing_random.raw_res', skipinitialspace=True)
	df_org['throughput'] = df_org['insertions'] / 1000000 / df_org['duration']
	del df_org['insertions']
	del df_org['duration']
	del df_org['dataset']
	qmax_df_org = df_org[(df_org['type'] == 'AmortizedQMax')]
	las_qmax_df_org = df_org[(df_org['type'] == 'AmortizedSampledQMaxV2')]
	las9_qmax_df_org = df_org[(df_org['type'] == 'AmortizedSampledQMaxV2U1')]
	ylabel = 'CPU Throughput [MPPS]'
	y_max=400
	x_min= 0.001
	x_max = 4
	save_dir="../figures/"
	legend_loc='upper left'
	legend_anchor=(0,1)
else:
	df_org = pd.read_csv('../ovs_results/graph3_' + G + '.txt', skipinitialspace=True)
	ylabel = 'OVS Throughput [Gbps]'
	df_org['throughput'] = df_org['throughput'] / 1000000000
	qmax_df_org = df_org[(df_org['type'] == 'QMax')]
	las_qmax_df_org = df_org[(df_org['type'] == 'Sampled QMax V2')]
	las9_qmax_df_org = df_org[(df_org['type'] == 'Sampled QMax')]
	y_max=int(G[:-1]) + 1
	x_min = 0.01
	x_max = 1
	save_dir="../ovs_figures/"
	legend_loc='lower right'
	legend_anchor=(1,0)

x_mid = 2**((math.log(x_max,2) + math.log(x_min,2))/2)

heap_df_org = df_org[(df_org['type'] == 'Heap')]
skiplist_df_org = df_org[(df_org['type'] == 'SkipList')]
ovs_df_org = df_org[(df_org['type'] == 'OVS')]
del ovs_df_org['size']

sizes=map(lambda s : 10**int(s), sys.argv[2:])
for s, oc, hc, sc, qc, las_qc, las9_qc in itertools.izip(sizes, ocs, hcs, scs, qcs,las_qcs, las9_qcs):
	fig, ax1 = plt.subplots(1,1)
#	ax1.grid()
	ax1.set_ylabel(ylabel, fontsize=20)
	ax1.set_xlabel(r'Performance parameter ($\gamma$)', fontsize=20)

	#q-MAX prepration
        qmax_size_df = qmax_df_org[(qmax_df_org['size'] == s)]
        
        del qmax_size_df['size']
        gq = qmax_size_df.groupby(['gamma']).agg(['mean', 'std'])
        gq = gq.reset_index()
	sg=len(gq['gamma'])
	
        las_qmax_size_df = las_qmax_df_org[(las_qmax_df_org['size'] == s)]
        del las_qmax_size_df['size']
        las_gq = las_qmax_size_df.groupby(['gamma']).agg(['mean', 'std'])
        las_gq = las_gq.reset_index()
        
        las9_qmax_size_df = las9_qmax_df_org[(las9_qmax_df_org['size'] == s)]
        del las9_qmax_size_df['size']
        las9_gq = las9_qmax_size_df.groupby(['gamma']).agg(['mean', 'std'])
        las9_gq = las9_gq.reset_index()

	#OVS
        #ax1.axhline(y=ovs_df_org.mean()['throughput'], label='OVS', markersize=10, color=oc['color'], marker=oc['marker'], linestyle=oc['linestyle'])
	if G != "CPU":
		#ax1.errorbar(gq['gamma'],[ovs_df_org.mean()['throughput']]*sg, yerr=[ovs_df_org.std()['throughput']]*sg, label='OVS', fmt='o', linewidth=2, capsize=10, elinewidth=3, capthick=1.5 ,markersize=10, color=oc['color'], marker=oc['marker'])
		ax1.errorbar(x_mid,ovs_df_org.mean()['throughput'], yerr=ovs_df_org.std()['throughput'], label='OVS', fmt='o', linewidth=2, capsize=1000, elinewidth=3, capthick=1.5 ,markersize=10, color=oc['color'], marker=oc['marker'])


	# Heap
	heap_size_df = heap_df_org[(heap_df_org['size'] == s)]
	del heap_size_df['size']
	gh = heap_size_df.mean()
        #ax1.axhline(y=gh['throughput'], label='Heap', markersize=10, color=hc['color'], marker=hc['marker'], linestyle=hc['linestyle'])
        #ax1.errorbar(gq['gamma'],[heap_size_df.mean()['throughput']]*sg, yerr=[heap_size_df.std()['throughput']]*sg, label='Heap', fmt='o', linewidth=2, capsize=10, elinewidth=3, capthick=1.5, markersize=10, color=hc['color'], marker=hc['marker'])
	#ax1.errorbar(x_mid,heap_size_df.mean()['throughput'], yerr=heap_size_df.std()['throughput'], label='Heap', fmt='o', linewidth=2, capsize=1000, elinewidth=3, capthick=1.5, markersize=10, color=hc['color'], marker=hc['marker'])


	# Skiplist
	skiplist_size_df = skiplist_df_org[(skiplist_df_org['size'] == s)]
	del skiplist_size_df['size']
	gs = skiplist_size_df.mean()
        #ax1.axhline(y=gs['throughput'], label='SkipList', markersize=10, color=sc['color'], marker=sc['marker'], linestyle=sc['linestyle'])
        #ax1.errorbar(gq['gamma'],[skiplist_size_df.mean()['throughput']]*sg, yerr=[skiplist_size_df.std()['throughput']]*sg, label='SkipList', fmt='o', linewidth=2, capsize=10, elinewidth=3, capthick=1.5, markersize=10, color=sc['color'], marker=sc['marker'])
	#ax1.errorbar(x_mid,skiplist_size_df.mean()['throughput'], yerr=skiplist_size_df.std()['throughput'], label='SkipList', fmt='o', linewidth=2, capsize=1000, elinewidth=3, capthick=1.5, markersize=10, color=sc['color'], marker=sc['marker'])

	#q-MAX
        #ax1.plot(gq.mean(), label='q-MAX', markersize=10)
	ax1.errorbar(gq['gamma'],gq['throughput']['mean'], yerr=gq['throughput']['std'], label='q-MAX',  fmt='o', linewidth=2, capsize=10, elinewidth=3, capthick=1.5, markersize=10, color=qc['color'], marker=qc['marker'])
	
	#las q-Max 
	ax1.errorbar(las_gq['gamma'],las_gq['throughput']['mean'], yerr=las_gq['throughput']['std'], label='Sampled q-MAX Las-Vegas V2',  fmt='o', linewidth=2, capsize=5, elinewidth=3, capthick=1.5, markersize=10, color=las_qc['color'], marker=las_qc['marker'])
	
	ax1.errorbar(las9_gq['gamma'],las9_gq['throughput']['mean'], yerr=las9_gq['throughput']['std'], label='Sampled q-MAX Las-Vegas MySuggestion',  fmt='o', linewidth=2, capsize=5, elinewidth=3, capthick=1.5, markersize=10, color=las9_qc['color'], marker=las9_qc['marker'])


	#Plot configurations
        handles, labels = ax1.get_legend_handles_labels()
	new_handles = []
	for h in handles:
		#only need to edit the errorbar legend entries
		if isinstance(h, container.ErrorbarContainer):
			new_handles.append(h[0])
		else:
			new_handles.append(h)
        ax1.legend(new_handles, labels, prop={'size': 16}, loc=legend_loc, bbox_to_anchor=legend_anchor, ncol=1)
        fig.tight_layout()

	ax1.set_xticks(np.arange(sg))
	plt.tick_params(axis='both', which='major', labelsize=16)
	plt.xlim(xmin=x_min, xmax=x_max)
	plt.ylim(ymin=0, ymax=y_max)
	ax1.set_xscale('log', basex=2)
	for item in ([ax1.xaxis.label, ax1.yaxis.label] +
	             ax1.get_xticklabels() + ax1.get_yticklabels()):
	    item.set_fontsize(24)
	#plt.show()
	plt.savefig(save_dir+'graph3-'+G+'-q='+str([s])+'.pdf', bbox_inches='tight')
