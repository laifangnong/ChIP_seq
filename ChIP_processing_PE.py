import sys,os
#import difflib
from  LCS_name_R_1_2 import longest_common_subsequence
from subprocess import run
from subprocess import call
from itertools import product
import subprocess
import time
class ChIP_seq_analyses:
	def __init__(self,fastq_dir):
		self.fastq_dir=fastq_dir
		self.all_pair=[]
		self.outputdir1='/'.join(self.fastq_dir.split('/')[0:-2])+'/step1_fastqc'
		self.outputdir2='/'.join(self.fastq_dir.split('/')[0:-2])+'/step2_trim_galore'
		self.outputdir3='/'.join(self.fastq_dir.split('/')[0:-2])+'/step3_bowtie2'
		self.logdir3='/'.join(self.fastq_dir.split('/')[0:-2])+'/log3_bowtie2'
		self.outputdir4='/'.join(self.fastq_dir.split('/')[0:-2])+'/step4_rm_cut'
		self.outputdir5='/'.join(self.fastq_dir.split('/')[0:-2])+'/step5_picard_polyclonal'
		self.outputdir6='/'.join(self.fastq_dir.split('/')[0:-2])+'/step6_rm_poly_random_chrM'
		self.outputdir7='/'.join(self.fastq_dir.split('/')[0:-2])+'/step7_sam2bam'
		self.outputdir8='/'.join(self.fastq_dir.split('/')[0:-2])+'/step8_bamCoverage'
		self.outputdir9='/'.join(self.fastq_dir.split('/')[0:-2])+'/step9_header'
		
	def step1_fastqc(self):
		fastq_dir=os.listdir(self.fastq_dir)
		
		if not os.path.exists(self.outputdir1):
			os.makedirs(self.outputdir1)
		else:
			pass
		range_fastq_dir=range(len(fastq_dir))
		procs1=[]
		for i in range_fastq_dir:
			for j in range_fastq_dir:
				if longest_common_subsequence(fastq_dir[i],fastq_dir[j])==1:
					if ('1.') in fastq_dir[i] and ('2.') in fastq_dir[j]:
						if not sorted([fastq_dir[i],fastq_dir[j]]) in self.all_pair:
							self.all_pair.append(sorted([fastq_dir[i],fastq_dir[j]]))
		for x in range(len(self.all_pair)):
			#step1_fastqc
			cmd1='nohup fastqc -o '+self.outputdir1+' '+self.fastq_dir+'/'+self.all_pair[x][0]+' '+self.fastq_dir+'/'+self.all_pair[x][1]+' >fastqc_run.log &'
			proc=subprocess.Popen(cmd1,shell=True,stdout=subprocess.PIPE)
			procs1.append(proc)
		for x in procs1:
			x.communicate()
			
			
			
	def step2_trim_galore(self):
		procs2=[]
		print (self.outputdir2)
		if not os.path.exists(self.outputdir2):
			os.makedirs(self.outputdir2)
		else:
			pass	
		for x in range(len(self.all_pair)):
			#step2_trim_galore
			cmd2='nohup trim_galore -q 25 --stringency 5 --fastqc --phred33 -dont_gzip --retain_unpaired -r1 35 -r2 35 --length 34 -o '+self.outputdir2+' --paired '+self.fastq_dir+'/'+self.all_pair[x][0]+' '+self.fastq_dir+'/'+self.all_pair[x][1]+' &'		
			proc=subprocess.Popen(cmd2,shell=True,stdout=subprocess.PIPE)
			procs2.append(proc)
		for x in procs2:
			x.communicate()
			
	def step3_bowtie2_alignment(self):		
		if not os.path.exists(self.outputdir3):
			os.makedirs(self.outputdir3)
		else:
			pass
		if not os.path.exists(self.logdir3):
			os.makedirs(self.logdir3)
		else:
			pass
		procs3=[]
		step2_outputdir=os.listdir(self.outputdir2)
		range_step2_outputdir=range(len(step2_outputdir))
		all_fastq_step2=[]
		all_pair_step2=[]
		all_prefix_step2=[]
		for i in range_step2_outputdir:
			if 'val' in step2_outputdir[i] and 'unpaired' not in step2_outputdir[i]:
				if step2_outputdir[i].endswith('fq') or step2_outputdir[i].endswith('fastq'):
					all_fastq_step2.append(step2_outputdir[i])
		range_step2_output=range(len(all_fastq_step2))
		for x,y in product(all_fastq_step2,all_fastq_step2):
			if longest_common_subsequence(x,y)==2 and sorted([x,y]) not in all_pair_step2:
				if '_r1_' and '_1' in x and '_r2_' and '_2' in y:
					all_pair_step2.append(sorted([x,y]))
					all_prefix_step2.append (sorted([x,y])[0].split('r1')[0][0:-1])
					cmd3='nohup bowtie2 -t -q -p 4 -N 1 -L 25 -X 2000 --no-mixed --no-discordant -x /data3/fangnong/genome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome -1 '+self.outputdir2+'/'+sorted([x,y])[0]+' -2 '+self.outputdir2+'/'+sorted([x,y])[1]+' > '+self.outputdir3+'/'+ sorted([x,y])[0].split('r1')[0][0:-1]+'.sam'+' 2>>'+self.logdir3+'/'+sorted([x,y])[0].split('r1')[0][0:-1]+'_mapping.log '
					proc=subprocess.Popen(cmd3,shell=True,stdout=subprocess.PIPE)
					procs3.append(proc)
		for x in procs3:
			x.communicate()
		
		
		

		
	def step4_rm_cut(self):	
		if not os.path.exists(self.outputdir4):
			os.makedirs(self.outputdir4)
		else:
			pass
		
		bowtie2_out=os.listdir(self.outputdir3)
		range_bowtie2_out=range(len(bowtie2_out))
		procs4=[]
		for i in range_bowtie2_out:
			cmd4='perl /data/jingyi/scripts/a_batch/rm_cut.pl '+self.outputdir3+'/'+bowtie2_out[i]+' > '+self.outputdir4+'/'+bowtie2_out[i][0:-4]+'.cleaned.sam'
			print (cmd4)
			proc=subprocess.Popen(cmd4,shell=True,stdout=subprocess.PIPE)
			procs4.append(proc)
		for x in procs4:
			x.communicate()

	def step5_picard_rmdup(self):
		if not os.path.exists(self.outputdir5):
			os.makedirs(self.outputdir5)
		else:
			pass
		procs5=[]
		rm_cut_out=os.listdir(self.outputdir4)
		range_rm_cut_out=range(len(rm_cut_out))
		for i in range_rm_cut_out:
			cmd5='bash /data4/bofeng/scripts/a_batch/rm_polyclonal.sh /data/jingyi/sharedata/genome/mm9/mm9 ' +self.outputdir4+'/'+rm_cut_out[i][0:-4]
			print (cmd5)
			proc=subprocess.Popen(cmd5,shell=True,stdout=subprocess.PIPE)
			procs5.append(proc)
		for x in procs5:
			x.communicate()
		procs6=[]
		rm_cut_out=os.listdir(self.outputdir4)
		range_rm_cut_out=range(len(rm_cut_out))
		for i in range_rm_cut_out:
			if rm_cut_out[i].endswith('.bam'):
				cmd6='mv '+self.outputdir4+'/'+rm_cut_out[i]+' '+self.outputdir5+'/'+rm_cut_out[i]
				proc=subprocess.Popen(cmd6,shell=True,stdout=subprocess.PIPE)
				procs6.append(proc)
		for x in procs6:
			x.communicate()
			
	def step6_filter_dup_random_chrM(self):
		if not os.path.exists(self.outputdir6):
			os.makedirs(self.outputdir6)
		else:
			pass
		procs7=[]
		rm_poly_out=os.listdir(self.outputdir5)
		range_rm_poly_out=range(len(rm_poly_out))
		for i in range_rm_poly_out:
			if rm_poly_out[i].endswith('cleaned.bam'):
				cmd7='samtools view -h '+self.outputdir5+'/'+rm_poly_out[i]+' > '+self.outputdir6+'/'+rm_poly_out[i][0:-4]+'.sam'
				proc=subprocess.Popen(cmd7,shell=True,stdout=subprocess.PIPE)
				procs7.append(proc)
		for x in procs7:
			x.communicate()
			
		procs8=[]
		rm_poly_out=os.listdir(self.outputdir6)
		range_rm_poly_out=range(len(rm_poly_out))
		for i in range_rm_poly_out:
			if rm_poly_out[i].endswith('.sam'):
				cmd8='python /data3/fangnong/script/ChIP_sam_filter.py '+self.outputdir6+'/'+rm_poly_out[i]
				proc=subprocess.Popen(cmd8,shell=True,stdout=subprocess.PIPE)
				procs8.append(proc)
		for x in procs8:
			x.communicate()
		
	
	
	
	def step7_sam2bam_index(self):
		if not os.path.exists(self.outputdir7):
			os.makedirs(self.outputdir7)
		else:
			pass
		procs9=[]
	#	procs10=[]
		filter_out=os.listdir(self.outputdir6)
		range_filter_out=range(len(filter_out))
		for i in range_filter_out:
			if filter_out[i].endswith('pmu.sam'):
				cmd9='samtools view -b -S '+self.outputdir6+'/'+filter_out[i]+' > '+self.outputdir7+'/'+filter_out[i][0:-4]+'.bam'
				print(cmd9)
				proc=subprocess.Popen(cmd9,shell=True,stdout=subprocess.PIPE)
				procs9.append(proc)
			#	cmd10='samtools index '+self.outputdir7+'/'+filter_out[i][0:-4]+'.bam'
			#	proc=subprocess.Popen(cmd10,shell=True,stdout=subprocess.PIPE)
			#	procs10.append(proc)
		for x in procs9:
			x.communicate()
		#for x in procs10:
		#	x.communicate()
	def step7_1_index(self):
		procs10=[]
		filter_out=os.listdir(self.outputdir6)			
		range_filter_out=range(len(filter_out))
		for i in range_filter_out:
			if filter_out[i].endswith('pmu.sam'):
				cmd10='samtools index '+self.outputdir7+'/'+filter_out[i][0:-4]+'.bam'
				proc=subprocess.Popen(cmd10,shell=True,stdout=subprocess.PIPE)
				procs10.append(proc)
		for x in procs10:
			x.communicate()
	#def step8_bam_index(self):
	
	def step8_bamCoverage_bw(self):
		if not os.path.exists(self.outputdir8):
			os.makedirs(self.outputdir8)
		else:
			pass
		procs11=[]
		filter_out=os.listdir(self.outputdir7)
		range_bam_out=range(len(filter_out))
		for i in range_bam_out:
			if filter_out[i].endswith('.bam'):
				cmd11='bamCoverage -b '+self.outputdir7+'/'+filter_out[i]+' -o '+self.outputdir8+'/'+filter_out[i][0:-4]+'.bw'+' -bs 100 --normalizeUsing RPKM'
				proc=subprocess.Popen(cmd11,shell=True,stdout=subprocess.PIPE)
				procs11.append(proc)
		for x in procs11:
			x.communicate()
	
	def step9_mv_webfolder(self):
		webfolder='/webhtml/fangnong/'+time.strftime("%Y_%m_%d",time.localtime()) 
		if not os.path.exists(webfolder):
			os.makedirs(webfolder)
		else:
			pass
		procs12=[]
		bw_out=os.listdir(self.outputdir8)
		range_bw_out=range(len(bw_out))
		for i in range_bw_out:
			if bw_out[i].endswith('.bw'):
				cmd12='cp '+self.outputdir8+'/'+bw_out[i]+' '+webfolder+'/'+bw_out[i]
				proc=subprocess.Popen(cmd12,shell=True,stdout=subprocess.PIPE)
				procs12.append(proc)
		for x in procs12:
			x.communicate()
				
	
	def step10_header(self):
		if not os.path.exists(self.outputdir9):
			os.makedirs(self.outputdir9)
		else:
			pass
		header=open(self.outputdir9+'/header.txt','w')
		bw_out=os.listdir(self.outputdir8)
		range_bw_out=range(len(bw_out))
		for i in range_bw_out:
			if bw_out[i].endswith('.bw'):
				header.write('track type=bigWig name='+bw_out[i][0:-3]+' description='+bw_out[i][0:-3]+' bigDataUrl=http://166.111.156.41/fangnong/'+time.strftime("%Y_%m_%d",time.localtime())+'/'+bw_out[i]+' viewLimits=3:12 visibility=2 windowingFunction=maximum maxHeightPixels=30 autoScale=off color=0,0,123'+'\n')
				
	#def step12_back_up_rm(self):
	
		#import time
		#webfolder='/webhtml/fangnong/'+time.strftime("%Y_%m_%d",time.localtime()) 
		#if not os.path.exists(webfolder):
		#	os.makedirs(webfolder)
		#else:
		#	pass
		#cmd_cd='cd '+webfolder
		#run (cmd_cd,shell=True)
		#procs4=[]
		#for i in range_bowtie2_out:
		#	cmd5='nohup bash sam2bw_chip_total.sh mm9 '+outputdir3+'/'+bowtie2_out[i]+' '+webfolder+'&'
		#	proc=subprocess.Popen(cmd5,shell=True,stdout=subprocess.PIPE)
		#	procs4.append(proc)
		#	print (cmd5)
		#for x in procs4:
		#	x.communicate()

C=ChIP_seq_analyses(sys.argv[1])	
#C.step1_fastqc()
#C.step2_trim_galore()
#C.step3_bowtie2_alignment()
#C.step4_rm_cut()
#C.step5_picard_rmdup()
#C.step6_filter_dup_random_chrM()
#C.step7_sam2bam_index()
#C.step7_1_index()
C.step8_bamCoverage_bw()
C.step9_mv_webfolder()
C.step10_header()
			

