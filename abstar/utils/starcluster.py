#!/usr/bin/env python
# filename: starcluster.py

#
# Copyright (c) 2015 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


from __future__ import print_function

import sys
import time
import logging
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
    raid = '10' if len(drives) == 24 else '0'
    log.write('\nBuilding a RAID{} array using {} instance storage drives on master node...'.format(raid, len(drives)))
    log.flush()
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
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    volume_size = get_volume_size(data_dir)
    log.write('Done.\nArray contains {} drives, mounted at {}\n'.format(len(drives), data_dir))
    if volume_size:
        log.write('{}B of space is available.\n'.format(volume_size))


def get_volume_size(volume):
    cmd = 'df -h'
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    for line in stdout.split('\n')[1:]:
        device, size, used, available, percent, mountpoint = line.split()
        if mountpoint == volume:
            return available
    return None


def share_data_array(data_dir, nodes, log):
    node_names = [n[0] for n in nodes]
    add_nodes_to_export(data_dir, node_names, log)
    mount_array_on_nodes(data_dir, node_names, log)


def add_nodes_to_export(data_dir, nodes, log):
    log.write('Adding nodes to /etc/exports on master node...\n')
    for node in nodes:
        cmd = 'echo "{} {}(async,no_root_squash,no_subtree_check,rw)" | sudo tee -a /etc/exports'.format(data_dir, node)
        p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
    p = sp.Popen('sudo /etc/init.d/nfs-kernel-server start', shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()


def mount_array_on_nodes(data_dir, nodes, log):
    log.write('Mounting NFS share on each node:\n')
    cmd = "sudo mkdir {0}\
           && sudo mount master:{0} {0}\
           && sudo chmod 777 {0}".format(data_dir)
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


def launch_celery_workers(nodes, log):
    node_names = [n[0] for n in nodes if n[0].lower() != 'master']
    cmd = 'cd /abstar\
           && chmod 777 -R /abstar\
           && sudo -u sgeadmin celery -A utils.queue.celery worker -l info --detach'
    for node in node_names:
        try:
            log.write(node + '\n')
            run_remote_cmd(node, cmd)
        except:
            log.write('FAILED\n')
            sys.exit()



def print_start_info(log):
    log.write('\n')
    log.write('\n')
    log.write('========================================\n')
    log.write('Configuring StarCluster\n')
    log.write('========================================\n')
    log.write('\n')



def run(data_dir='/data', log=sys.stdout):
    print_start_info(log)
    log.write('Identifying all cluster nodes...')
    nodes = identify_nodes()
    log.write('{} nodes identified\n'.format(len(nodes)))
    logging.info('Cluster size: {} nodes'.format(len(nodes)))

    build_data_array(data_dir, log)

    log.write('\nSharing data array with all nodes via NFS...\n')
    share_data_array(data_dir, nodes, log)
    logging.info('NFS sharing configured for all nodes.')

    log.write('\nLaunching Celery workers on all worker nodes...\n')
    launch_celery_workers(nodes, log)
    logging.info('Celery workers started all nodes.')

    log.write('\nDone. StarCluster is now properly configured for an AbStar run.\n\n')



if __name__ == '__main__':
    # standalone args
    parser = argparse.ArgumentParser("")
    parser.add_argument('-d', '--data', dest='data_dir', default='/data', help="The location in which to mount the RAID data array. Default is '/data'.")
    parser.add_argument('-l', '--log', dest='log', default=None, help="The log file. If not provided, log info will be written to sys.stdout.")
    args = parser.parse_args()

    args.log = args.log if args.log else sys.stdout

    run(data_dir=args.data_dir, log=args.log)
