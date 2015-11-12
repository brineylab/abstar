#!/usr/bin/env python
# filename: starcluster.py

from __future__ import print_function

import os
import sys
import time
import string
import argparse
import subprocess as sp

import paramiko


def identify_nodes():
	nodes = []
	with open('/etc/hosts', 'r') as f:
		lines = iter(f.readlines())
		line = lines.next()
		while 'master' not in line:
			line = lines.next()
		try:
			while True:
				ip, name = line.rstrip('\n').split()
				nodes.append((name, ip))
				line = lines.next()
		except StopIteration:
			pass
	return nodes


def identify_drives():
	cmd = 'ls /dev/'
	p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	stdout, stderr = p.communicate()
	drives = ['/dev/{}'.format(s) for s in stdout.split() if s.startswith('xvda') and s != 'xvda']
	return drives


def build_data_array(data_dir, log):
	drives = identify_drives()
	raid = '10' if len(drives) >= 4 and len(drives) % 2 == 0 else '0'
	log.write('\nBuilding a RAID{} array using {} instance storage drives on master node...\n'.format(raid, len(drives)))
	cmd = "sudo mdadm --verbose --create /dev/md0 --level={3} --chunk=256 --raid-devices={0} {1}\
		   && sudo dd if=/dev/zero of=/dev/md0 bs=512 count=1\
		   && sudo pvcreate /dev/md0\
		   && sudo vgcreate vg0 /dev/md0\
		   && sudo lvcreate -l 100%vg -n data vg0\
		   && sudo mke2fs -t ext4 -F /dev/vg0/data\
		   && sudo mkdir {2}\
		   && echo '/dev/vg0/data {2} ext4 defaults,auto,noatime,noexec 0 0' | sudo tee -a /etc/fstab\
		   && sudo mount {2}\
		   && sudo chmod 777 {2}\
		   ".format(len(drives), ' '.join(drives), data_dir, raid)
	os.system(cmd)
	log.write('Array contains {} drives, mounted at {}\n'.format(len(drives), data_dir))


def share_data_array(data_dir, nodes, log):
	node_names = [n[0] for n in nodes]
	add_nodes_to_export(data_dir, node_names, log)
	mount_array_on_nodes(data_dir, node_names, log)


def add_nodes_to_export(data_dir, nodes, log):
	log.write('Adding nodes to /etc/exports on master node...\n')
	for node in nodes:
		os.system('echo "\n{} {}(async,no_root_squash,no_subtree_check,rw)" | sudo tee -a /etc/exports'.format(data_dir, node))
	os.system('sudo /etc/init.d/nfs-kernel-server start')


def mount_array_on_nodes(data_dir, nodes, log):
	log.write('Mounting NFS share on each node:\n')
	cmd = "sudo mkdir {0}\
		   && sudo chmod 777 {0}\
		   && sudo mount master:{0} {0}".format(data_dir)
	for node in nodes:
		log.write(node + '\n')
		run_remote_cmd(node, cmd)


def run_remote_cmd(node, cmd):
		'''
		Remotely runs 'cmd' on 'node' and blocks until 'cmd' completes.
		Returns 'cmd' stdout.
		'''
		ssh = paramiko.SSHClient()
		ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		ssh.connect(node)
		stdin, stdout, stderr = ssh.exec_command(cmd)
		channel = stdout.channel
		while not channel.exit_status_ready():
			time.sleep(1)
		return stdout.read()


def build_rabbitmq_cluster(nodes, log):
	# node_names = [n[0] for n in nodes if n[0].lower() != 'master']
	master_ip = nodes[0][1]
	node_ips = [n[1]for n in nodes if n[0].lower() != 'master']
	cmd = 'sudo rabbitmqctl stop_app \
		  && sudo rabbitmqctl join_cluster rabbit@{} \
		  && sudo rabbitmqctl start_app'.format(master_ip)
	for node in node_ips:
		try:
			log.write(node + '\n')
			stdout = run_remote_cmd(node, cmd)
			print(stdout)
		except:
			log.write('FAILED\n')
			sys.exit()


def launch_celery_workers(nodes, log):
	node_names = [n[0] for n in nodes]
	# node_ips = [n[1]for n in nodes if n[0].lower() != 'master']
	cmd = 'cd /celery/celery_tests\
		  && chmod 777 -R /celery/celery_tests/\
		  && sudo -u sgeadmin celery -A utils.celery.celery worker -l info --detach'
	for node in node_names:
		try:
			log.write(node + '\n')
			run_remote_cmd(node, cmd)
		except:
			log.write('FAILED\n')
			sys.exit()



def run(data_dir='/data', log=sys.stdout):
	log.write('\nIdentifying all cluster nodes...\n')
	nodes = identify_nodes()

	# build_data_array(data_dir, log)

	# log.write('\nSharing data array with all nodes via NFS...\n')
	# share_data_array(data_dir, nodes, log)

	# log.write('\nJoining all nodes into a RabbitMQ cluster...\n')
	# build_rabbitmq_cluster(nodes, log)

	log.write('\nLaunching Celery workers on all nodes...\n')
	launch_celery_workers(nodes, log)

	log.write('\nDone. StarCluster is now properly configured for an AbAnalyze run.\n\n')



if __name__ == '__main__':
	# standalone args
	parser = argparse.ArgumentParser("")
	parser.add_argument('-d', '--data', dest='data_dir', default='/data', help="The location in which to mount the RAID data array. Default is '/data'.")
	parser.add_argument('-l', '--log', dest='log', default=None, help="The log file. If not provided, log info will be written to sys.stdout.")
	args = parser.parse_args()

	args.log = args.log if args.log else sys.stdout
	
	run(data_dir=args.data_dir, log=args.log)
