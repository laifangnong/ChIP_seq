import subprocess
def wget_ref(ref_SRR_list):
	ref_SRR=open(ref_SRR_list,'r')
	procs=[]
	for i in ref_SRR:
		if i.startswith('SRR'):
			cmd='nohup wget -P /data3/fangnong/ref_public_ChIP/sra  -r ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/'+i.split('\t')[0][0:6]+'/'+i.split('\t')[0]+'/ &'
			print (cmd)
			proc=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
			procs.append(proc)
	for x in procs:
		x.communicate()
import sys,os		
wget_ref(sys.argv[1])
from subprocess import run
def mv_to_root(rootdir):
	root_dir=os.listdir(rootdir)
	range_root_dir=range(len(root_dir))
	for i in range_root_dir:
		if os.path.isdir(rootdir+'/'+root_dir[i]):
			mv_to_root(rootdir+'/'+root_dir[i])
		else:
			if os.path.isfile(rootdir+'/'+root_dir[i]):
				cmd='mv '+rootdir+'/'+root_dir[i]+' /data3/fangnong/ref_public_ChIP/sra_new'
				run(cmd,shell=True)
mv_to_root(sys.argv[2])
#python wget_ref.py /data3/fangnong/ref_public_ChIP/list.txt /data3/fangnong/ref_public_ChIP/sra
