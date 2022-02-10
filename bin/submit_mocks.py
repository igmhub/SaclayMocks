#!/usr/bin/env python
from __future__ import division, print_function
import os, sys
import argparse
import subprocess
import healpy as hp
import numpy as np
from SaclayMocks import util


"""
This code is highly inspired from picca_nersc_submit.py (https://github.com/igmhub/picca)
This code produces all the need bash and sbatch script to run a full production.
There are several parameters to provide with argparse. Some other can be tunned
manualy in the main function.

At the end, this code produces a bash file: submit.sh, for each mock realisation.
You need to run this bash file in order to send the jobs to cori nodes.
08/04/2019 - Thomas Etourneau
thomas.etouneau@cea.fr
"""


def get_header(mock_args, sbatch_args, jobname):
    '''
    returns the standard header for running the analyses.
    Input:
        - time (string): requested maximum runtime, format: hh:mm:ss
        - name (string): name of the job
        - email (string, optional): email address to send updates
        - queue (string, optional): name of the submission queue
    Returns:
        a string holding the header
    Function from picca_nersc_submit.py,
    https://github.com/igmhub/picca
    '''
    time = sbatch_args['time_{}'.format(jobname)]
    node = sbatch_args['nodes_{}'.format(jobname)]
    name = sbatch_args['name_{}'.format(jobname)]
    email = sbatch_args['email']
    queue = sbatch_args['queue_{}'.format(jobname)]
    account = sbatch_args['account']
    header = ""
    header += "#!/bin/bash\n"
    header += "#SBATCH -N {}\n".format(node)
    header += "#SBATCH -C haswell\n"
    header += "#SBATCH -q {}\n".format(queue)
    header += "#SBATCH -J {}\n".format(name)
    if email != None:
        header += "#SBATCH --mail-user={}\n".format(email)
        header += "#SBATCH --mail-type=ALL\n"
    header += "#SBATCH -t {}\n".format(time)
    header += "#SBATCH -L project\n"
    header += "#SBATCH -A {}\n".format(account)
    if mock_args['burst_buffer'] and jobname != 'pk':
        header += "#DW persistentdw name={}\n".format(mock_args['bb_name'])
    header += "#OpenMP settings:\n"
    header += "export OMP_NUM_THREADS=1\n"
    header += "start=$SECONDS\n"
    return header


def create_reservation(mock_args):
    '''
    This function creates a script that creates a Burst Buffer persistent reservation
    (To optimize the I/O on cori nodes)

    '''
    script = ""
    script += "#!/bin/bash\n"
    script += "#SBATCH -N 1\n"
    script += "#SBATCH -C haswell\n"
    script += "#SBATCH -J create_saclay_{}\n".format(mock_args['imock'])
    script += "#SBATCH -q regular\n"
    script += "#SBATCH -t 00:05:00\n"
    script += "#BB create_persistent name={name} capacity={size} access_mode=striped type=scratch\n".format(name=mock_args['bb_name'], size=mock_args['bb_size'])
    # script += "#DW persistentdw name={}\n".format(mock_args['bb_name'])
    # script += "#DW stage_in source={path} destination=$DW_PERSISTENT_STRIPED_{name}/pk type=directory\n".format(path=mock_args['dir_pk'], name=mock_args['bb_name'])
    # script += "#DW stage_in source={path} destination=$DW_PERSISTENT_STRIPED_{name}/mock_{i} type=directory\n".format(path=mock_args['base_dir'], name=mock_args['bb_name'], i=mock_args['imock'])

    # script += "echo $DW_PERSISTENT_STRIPED_{}\n".format(mock_args['bb_name'])
    # script += "echo '--'\n"
    # script += "ls -ahls $DW_PERSISTENT_STRIPED_{}\n".format(mock_args['bb_name'])
    # script += "echo '--'\n"
    # script += "ls -ls $DW_PERSISTENT_STRIPED_{}/pk\n".format(mock_args['bb_name'])
    # script += "echo '--'\n"
    # script += "ls -ls $DW_PERSISTENT_STRIPED_{}/mock_*\n".format(mock_args['bb_name'])
    # script += "echo '--'\n"
    # script += "echo 'number of files:'\n"
    # script += "ls -ls $DW_PERSISTENT_STRIPED_{}/* | wc -l \n".format(mock_args['bb_name'])
    script += "echo 'END'\n"

    filename = mock_args['run_dir']+'/create_reservation.sh'
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()


def stage_in(mock_args):
    '''
    This functions stages in the mock files on the burst buffer nodes.
    '''
    pkdir = 'dir_pk'
    if mock_args['burst_buffer']:
        pkdir += '_no_bb'
    script = "#!/bin/bash\n"
    script += "#SBATCH -N 1\n"
    script += "#SBATCH -C haswell\n"
    script += "#SBATCH -J stagein_saclay_{}\n".format(mock_args['imock'])
    script += "#SBATCH -q regular\n"
    script += "#SBATCH -t 00:05:00\n"
    script += "#DW persistentdw name={}\n".format(mock_args['bb_name'])
    script += "#DW stage_in source={path} destination=$DW_PERSISTENT_STRIPED_{name}/pk type=directory\n".format(path=mock_args[pkdir], name=mock_args['bb_name'])
    script += "#DW stage_in source={path} destination=$DW_PERSISTENT_STRIPED_{name}/mock_{i} type=directory\n".format(path=mock_args['base_dir'], name=mock_args['bb_name'], i=mock_args['imock'])

    script += "echo 'moving runs directory out (avoid stagin it out)'\n"
    script += '''
if [ ! -d "$DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard" ]; then
    mkdir $DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard
fi
'''.format(name=mock_args['bb_name'], i=mock_args['imock'])
    script += "mv $DW_PERSISTENT_STRIPED_{name}/mock_{i}/output/runs $DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard/runs\n".format(name=mock_args['bb_name'], i=mock_args['imock'])
    script += "mv $DW_PERSISTENT_STRIPED_{name}/mock_{i}/output/logs $DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard/logs\n".format(name=mock_args['bb_name'], i=mock_args['imock'])
    script += "echo $DW_PERSISTENT_STRIPED_{}\n".format(mock_args['bb_name'])
    script += "echo '--'\n"
    script += "ls -ahls $DW_PERSISTENT_STRIPED_{}\n".format(mock_args['bb_name'])
    script += "echo '--'\n"
    script += "ls -ls $DW_PERSISTENT_STRIPED_{}/pk\n".format(mock_args['bb_name'])
    script += "echo '--'\n"
    script += "ls -ls $DW_PERSISTENT_STRIPED_{}/mock_*\n".format(mock_args['bb_name'])
    script += "echo '--'\n"
    script += "ls -ls $DW_PERSISTENT_STRIPED_{}/mock_*/output\n".format(mock_args['bb_name'])
    script += "echo '--'\n"
    script += "ls -ls $DW_PERSISTENT_STRIPED_{}/mock_*_discard\n".format(mock_args['bb_name'])
    script += "echo '--'\n"
    script += "echo 'number of files:'\n"
    script += "ls -ls $DW_PERSISTENT_STRIPED_{}/* | wc -l \n".format(mock_args['bb_name'])
    script += "echo 'END'\n"

    filename = mock_args['run_dir']+'/stage_in.sh'
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()


def stage_out_dir(mock_args):
    '''
    This function returns the bash lines to stage_out the specified directories
    '''
    script = ''
    # stage_out everything
    if 'all' in mock_args['stage_out_dir']:
        script += "#DW stage_out source=$DW_PERSISTENT_STRIPED_{name}/mock_{i} destination={path}/mock_{i} type=directory\n".format(name=mock_args['bb_name'], i=mock_args['imock'], path=mock_args['mock_dir'])
        return script

    # stage_out only some particular directories
    # output dir
    # if 'output' in mock_args['stage_out_dir']:
    script += "#DW stage_out source=$DW_PERSISTENT_STRIPED_{name}/mock_{i}/output destination={path}/mock_{i}/output type=directory\n".format(name=mock_args['bb_name'], i=mock_args['imock'], path=mock_args['mock_dir'])
    # boxes
    if 'boxes' in mock_args['stage_out_dir']:
        for j in mock_args['chunkid']:
            script += "#DW stage_out source=$DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/boxes destination={path}/mock_{i}/chunk_{j}/boxes type=directory\n".format(name=mock_args['bb_name'], i=mock_args['imock'], j=j, path=mock_args['mock_dir'])
    # qso
    if 'qso' in mock_args['stage_out_dir']:
        for j in mock_args['chunkid']:
            script += "#DW stage_out source=$DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/qso destination={path}/mock_{i}/chunk_{j}/qso type=directory\n".format(name=mock_args['bb_name'], i=mock_args['imock'], j=j, path=mock_args['mock_dir'])
    # randoms
    if 'randoms' in mock_args['stage_out_dir']:
        for j in mock_args['chunkid']:
            script += "#DW stage_out source=$DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/randoms destination={path}/mock_{i}/chunk_{j}/randoms type=directory\n".format(name=mock_args['bb_name'], i=mock_args['imock'], j=j, path=mock_args['mock_dir'])
    # spectra
    if 'spectra' in mock_args['stage_out_dir']:
        for j in mock_args['chunkid']:
            script += "#DW stage_out source=$DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/spectra destination={path}/mock_{i}/chunk_{j}/spectra type=directory\n".format(name=mock_args['bb_name'], i=mock_args['imock'], j=j, path=mock_args['mock_dir'])
    # spectra_merged
    if 'spectra_merged' in mock_args['stage_out_dir']:
        for j in mock_args['chunkid']:
            script += "#DW stage_out source=$DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/spectra_merged destination={path}/mock_{i}/chunk_{j}/spectra_merged type=directory\n".format(name=mock_args['bb_name'], i=mock_args['imock'], j=j, path=mock_args['mock_dir'])
    # dla
    if 'dla' in mock_args['stage_out_dir']:
        for j in mock_args['chunkid']:
            script += "#DW stage_out source=$DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/dla.fits destination={path}/mock_{i}/chunk_{j}/ type=file\n".format(name=mock_args['bb_name'], i=mock_args['imock'], j=j, path=mock_args['mock_dir'])
    # random dla
    if 'dla_randoms' in mock_args['stage_out_dir']:
        for j in mock_args['chunkid']:
            script += "#DW stage_out source=$DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/dla_randoms.fits destination={path}/mock_{i}/chunk_{j}/ type=file\n".format(name=mock_args['bb_name'], i=mock_args['imock'], j=j, path=mock_args['mock_dir'])

    return script


def discard_boxes(mock_args):
    script = ''
    script += "echo 'Boxes will be discarded. Only boxk.npy will be staged out.'\n"
    script += 'timer=$SECONDS\n'
    script += '''
if [ ! -d "$DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard" ]; then
    mkdir $DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard
fi
'''.format(name=mock_args['bb_name'], i=mock_args['imock'])
    # script += 'mkdir $DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard\n'.format(name=mock_args['bb_name'], i=mock_args['imock'])
    for j in mock_args['chunkid']:
        # script += 'mkdir $DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard/boxes_{j}\n'.format(name=mock_args['bb_name'], i=mock_args['imock'], j=j)
        script += 'mv $DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/boxes $DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard/boxes_{j}\n'.format(name=mock_args['bb_name'], i=mock_args['imock'], j=j)
        script += 'mkdir $DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/boxes\n'.format(name=mock_args['bb_name'], i=mock_args['imock'], j=j)
        script += 'mv $DW_PERSISTENT_STRIPED_{name}/mock_{i}_discard/boxes_{j}/boxk.npy $DW_PERSISTENT_STRIPED_{name}/mock_{i}/chunk_{j}/boxes\n'.format(name=mock_args['bb_name'], i=mock_args['imock'], j=j)
    script += '''echo "Boxes moved. Took $(( SECONDS - timer )) s"\n'''
    return script


def stage_out(mock_args):
    '''
    This functions stages out the mock outputs written on the burst buffer nodes.
    '''
    script = "#!/bin/bash\n"
    script += "#SBATCH -N 1\n"
    script += "#SBATCH -C haswell\n"
    script += "#SBATCH -J stageout_saclay_{}\n".format(mock_args['imock'])
    script += "#SBATCH -q regular\n"
    script += "#SBATCH -t 00:05:00\n"
    script += "#DW persistentdw name={}\n".format(mock_args['bb_name'])
    # add the #DW stage_out lines for each directory to stage out
    script += stage_out_dir(mock_args)
    script += "echo 'This will stage out:'\n"
    script += "echo '{}'\n".format(mock_args['stage_out_dir'])
    if mock_args['discard_boxes']:
        script += discard_boxes(mock_args)
    script += "echo 'END'\n"

    filename = mock_args['run_dir']+'/stage_out.sh'
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()


def delete_reservation(mock_args):
    '''
    This function delete the BB persistent reservation
    '''
    script = ""
    script += "#!/bin/bash\n"
    script += "#SBATCH -N 1\n"
    script += "#SBATCH -C haswell\n"
    script += "#SBATCH -J delete_saclay_{}\n".format(mock_args['imock'])
    script += "#SBATCH -q regular\n"  # don't need debug since it's just to free space in BB nodes at the total end
    script += "#SBATCH -t 00:05:00\n"
    script += "#BB destroy_persistent name={name}\n".format(name=mock_args['bb_name'])
    script += "echo 'END'\n"

    filename = mock_args['run_dir']+'/delete_reservation.sh'
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()


def change_permission():
    aa #prov


def get_errors(code, start, pids=None, threads_num=None):
    '''
    returns the bash lines to get potential errors in job runs
    '''
    errors = """
error=0
i=<start>
for p in $<pids>; do
    if wait $p; then
	echo "<code>-$i OK"
    else"""
    if threads_num is not None:
        errors += """
        if [ $((i % {})) -eq 0 ]; then
            echo "Error in <code>-$i but continue"
        else
	    echo "Error in <code>-$i"
	    let "error=$error+1"
        fi""".format(threads_num)
    else:
        errors += '''
	echo "Error in <code>-$i"
	let "error=$error+1"'''
    errors += """
    fi
    let "i=$i+1"
done
if [ $error -ne 0 ]; then
    echo -e "==> Error in $error <code> ...   Abort!"
    exit 1
fi
"""

    if pids is None:
        pids = 'pids'
    errors = errors.replace('<code>', str(code))
    errors = errors.replace('<start>', str(start))
    errors = errors.replace('<pids>', pids)
    return errors


def run_bash_script(codename, mock_args, sbatch_args):
    '''
    Produce the script to run bash script on cori node.
    This is used to run python scripts:
    draw_aso.py
    make_spectra.py
    merge_spectra.py
    Args:
    - scriptname : string that specifies the name of the bash script. run_draw_qso.sh for instance
    - codename : string that specifies the name of the python code. draw_qso for instance
    '''
    script = """echo -e "*** Running <codename> <nslice> times ***"\n"""
    script += """pids=""\n"""
    for node in range(sbatch_args['nodes_chunk']):
        if mock_args['sbatch']:
            script += "srun "
            if mock_args['verbosity'] is not None:
                script += mock_args['verbosity']
            script += " -N 1 -n 1 -c 64 "
        script += "bash {loc}/run_<codename>-{i_chunk}-{node}.sh ".format(
            loc=mock_args['run_dir_chunk-{}'.format(mock_args['i_chunk'])], i_chunk=mock_args['i_chunk'], node=node)
        script += ">> {logdir}/run_chunks-{ichunk}.log &\n".format(
            logdir=mock_args['logs_dir'], ichunk=mock_args['i_chunk'])
        script += """pids+=" $!"\n"""
    script += get_errors("run_"+codename, 0)
    script += """echo -e "==> run_<codename>.sh done. $(( SECONDS - start )) s"\n"""
    script = script.replace('<codename>', str(codename))
    script = script.replace('<nslice>', str(mock_args['nslice']))
    return script


def pk(mock_args, sbatch_args):
    '''
    write a .sh file to submit a job to produce the power spectrum
    '''
    script = get_header(mock_args, sbatch_args, "pk")
    # Run interpolate_pk.py :
    script += """echo -e "*** Running interpolate_pk {threads} times ***"\n""".format(threads=sbatch_args['threads_pk'])
    script += """pids=""\n"""
    for job in range(sbatch_args['threads_pk']):
        if mock_args['use_time']:
            script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
        script += "interpolate_pk.py -NX {nx} -NY {ny} -NZ {nz} -i {job} -N {nslice} -outDir {path} -pixel {pixel} ".format(
            nx=mock_args['nx'], ny=mock_args['ny'], nz=mock_args['nz'], job=job,
            nslice=sbatch_args['threads_pk'], path=mock_args['dir_pk'], pixel=mock_args['pixel_size'])
        script += "&> {logdir}/interpolate_pk-{job}.log &\n".format(logdir=mock_args['dir_pk_logs'], job=job)
        script += """pids+=" $!"\n"""
    script += get_errors("interpolate_pk.py", 0)
    script += """echo -e "==> interpolate_pk done. $(( SECONDS - start )) s"\n"""

    # Run merge_pk.py :
    script += """echo -e "*** Running merge_pk ***"\n"""
    if mock_args['use_time']:
        script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
    script += "merge_pk.py -NX {nx} -NY {ny} -NZ {nz} -N {nslice} -outDir {path} -inDir {path} ".format(
        nx=mock_args['nx'], ny=mock_args['ny'], nz=mock_args['nz'],
        nslice=sbatch_args['threads_pk'], path=mock_args['dir_pk'])
    script += "&> {logdir}/merge_pk.log".format(logdir=mock_args['dir_pk_logs'])
    script += """
if [ $? -ne 0 ]; then
    echo -e "/!\ Error in merge_pk ...  Abort!"
    exit 1
else
    echo -e "==> merge_pk done. $(( SECONDS - start )) s"
fi
echo "----- END -----"
"""

    filename = mock_args['dir_pk_run']+'/run_pk.sh'
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()


def boxes(mock_args, sbatch_args):
    '''
    write a .sh file to submit a job to produce the different boxes
    '''
    script = get_header(mock_args, sbatch_args, "boxes")
    script += """echo "Running run_boxes.sh"\n"""
    script += """echo "command: make_boxes.py -NX {nx} -NY {ny} -NZ {nz} -nHDU {nslice} -PkDir {path_pk} -outDir {path_boxes} -ncpu {threads} -pixel {pixel} -rsd {rsd} {seed} "\n""".format(nx=mock_args['nx'], ny=mock_args['ny'], nz=mock_args['nz'], nslice=mock_args['nslice'], path_pk=mock_args['dir_pk'], threads=sbatch_args['threads_boxes'], path_boxes=mock_args['dir_boxes-{}'.format(mock_args['i_chunk'])], pixel=mock_args['pixel_size'], rsd=mock_args['rsd'], seed=mock_args['seed'])
    if mock_args['use_time']:
        script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
    if mock_args['sbatch']:
        script += "srun "
        if mock_args['verbosity'] is not None:
            script += mock_args['verbosity']
        script += " -N 1 -n 1 -c 64 "
    script += "make_boxes.py -NX {nx} -NY {ny} -NZ {nz} -nHDU {nslice} -PkDir {path_pk} -outDir {path_boxes} -ncpu {threads} -pixel {pixel} -rsd {rsd} {seed} ".format(nx=mock_args['nx'], ny=mock_args['ny'], nz=mock_args['nz'], nslice=mock_args['nslice'], path_pk=mock_args['dir_pk'], threads=sbatch_args['threads_boxes'], path_boxes=mock_args['dir_boxes-{}'.format(mock_args['i_chunk'])], pixel=mock_args['pixel_size'], rsd=mock_args['rsd'], seed=mock_args['seed'])
    script += "&> {path}/make_boxes.log \n".format(path=mock_args['logs_dir_chunk-{}'.format(mock_args['i_chunk'])])
    script += """
if [ $? -ne 0 ]; then
    echo "==> Error in make_boxes ...   Abort!"
    exit 1
else
    echo -e "==> make_boxes done. $(( SECONDS - start )) s"
fi
echo "----- END -----"
"""

    filename = mock_args['run_dir']+'/run_boxes-{}.sh'.format(mock_args['i_chunk'])
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()


def chunk(todo, mock_args, sbatch_args):
    '''
    Write a .sh file to submit jobs that produce the different chunks
    '''
    script = get_header(mock_args, sbatch_args, "chunk")
    if "qso" in todo:
        script += run_bash_script('draw_qso', mock_args, sbatch_args)
    if "randoms" in todo:
        script += run_bash_script('randoms', mock_args, sbatch_args)
    if "make_spectra" in todo:
        script += run_bash_script('make_spectra', mock_args, sbatch_args)
    if "merge_spectra" in todo:
        script += run_bash_script('merge_spectra', mock_args, sbatch_args)
    script += """echo "----- END -----"\n"""

    filename = mock_args['run_dir']+'/run_chunk-{}.sh'.format(mock_args['i_chunk'])
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()


def mergechunks(todo, mock_args, sbatch_args):
    '''
    Write a .sh file to submit jobs that write the mock outputs in desi format.
    The python scripts involve here gather the output from the different chunks
    '''
    script = get_header(mock_args, sbatch_args, "mergechunks")
    if "merge_qso" in todo:
        script += """echo -e "*** Running merge_qso ***"\n"""
        if mock_args['use_time']:
            script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
        script += "merge_qso.py -inDir {inpath} -outDir {outpath} -nside {nside} -nest {nest} -zmin {zmin} -zmax {zmax} ".format(inpath=mock_args['base_dir'], outpath=mock_args['out_dir'], nside=mock_args['nside'], nest=mock_args['nest'], zmin=mock_args['zmin'], zmax=mock_args['zmax'])
        script += "&> {path}/merge_qso.log &\n".format(path=mock_args['logs_dir_mergechunks'])
        script += "pid_qso=$!\n"

    if "merge_randoms" in todo:
        script += """echo -e "*** Running merge_qso for randoms ***"\n"""
        if mock_args['use_time']:
            script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
        script += "merge_qso.py -inDir {inpath} -outDir {outpath} -nside {nside} -nest {nest} -zmin {zmin} -zmax {zmax} -random True ".format(inpath=mock_args['base_dir'], outpath=mock_args['out_dir'], nside=mock_args['nside'], nest=mock_args['nest'], zmin=mock_args['zmin'], zmax=mock_args['zmax'])
        script += "&> {path}/merge_randoms.log &\n".format(path=mock_args['logs_dir_mergechunks'])
        script += "pid_rand=$!\n"

    if "merge_qso" in todo:
        script += """
if wait $pid_qso; then
    echo "merge_qso OK"
else
    echo "Error in merge_qso"
    exit 1
fi
"""
    if "merge_randoms" in todo:
        script += """
if wait $pid_rand; then
    echo "merge_randoms OK"
else
    echo "Error in merge_randoms"
    exit 1
fi
"""
    if "merge_qso" in todo or "merge_randoms" in todo:
        script += """echo -e "==> QSO catalogs done. $(( SECONDS - start )) s"\n"""

    if "compute_dla" in todo:
        script += """echo -e "*** Producing DLA ***"\n"""
        script += "pids_dla=''\n"
        for cid in mock_args['chunkid']:
            if mock_args['use_time']:
                script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
            script += "dla_saclay.py --input_path {base}/chunk_{i}/spectra_merged/ --output_path {base}/chunk_{i} --cell_size {pixel} --nmin {nmin} --nmax {nmax} {seed} ".format(base=mock_args['base_dir'], i=cid, pixel=mock_args['pixel_size'], nmin=mock_args['nmin'], nmax=mock_args['nmax'], seed=mock_args['seed'])
            if mock_args['nhi_low_cut'] is not None and mock_args['nhi_high_cut'] is not None:
                script += " --nhi-low-cut {cut1} --nhi-high-cut {cut2} ".format(cut1=mock_args['nhi_low_cut'], cut2=mock_args['nhi_high_cut'])
            script += "&> {path}/dla-{i}.log &\n".format(path=mock_args['logs_dir_mergechunks'], i=cid)
            script += """pids_dla+=" $!"\n"""

    if "dla_randoms" in todo:
        script += """echo -e "*** Producing DLA randoms ***"\n"""
        script += "pids_rand=''\n"
        for cid in mock_args['chunkid']:
            if mock_args['use_time']:
                script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
            script += "dla_saclay.py --input_path {base}/chunk_{i}/spectra_merged/ --output_path {base}/chunk_{i} --cell_size {pixel} --nmin {nmin} --nmax {nmax} {seed} -random True ".format(base=mock_args['base_dir'], i=cid, pixel=mock_args['pixel_size'], nmin=mock_args['nmin'], nmax=mock_args['nmax'], seed=mock_args['seed'])
            if mock_args['nhi_low_cut'] is not None and mock_args['nhi_high_cut'] is not None:
                script += " --nhi-low-cut {cut1} --nhi-high-cut {cut2} ".format(cut1=mock_args['nhi_low_cut'], cut2=mock_args['nhi_high_cut'])
            script += "&> {path}/dla_rand-{i}.log &\n".format(path=mock_args['logs_dir_mergechunks'], i=cid)
            script += """pids_rand+=" $!"\n"""

    if "compute_dla" in todo:
        script += get_errors("dla", 1, pids="pids_dla")
    if "dla_randoms" in todo:
        script += get_errors("dla_rand", 1, pids="pids_rand")
    if "compute_dla" in todo or "dla_randoms" in todo:
        script += """echo -e "==> dla_saclay done. $(( SECONDS - start )) s"\n"""

    if "merge_dla" in todo:
        script += """echo -e "*** Merging DLA ***"\n"""
        if mock_args['use_time']:
            script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
        script += "merge_dla.py -indir {base} -outdir {outpath} ".format(base=mock_args['base_dir'], outpath=mock_args['out_dir'])
        if mock_args['nhi_low_cut'] is not None and mock_args['nhi_high_cut'] is not None:
            script += " --nhi-low-cut {cut1} --nhi-high-cut {cut2} ".format(cut1=mock_args['nhi_low_cut'], cut2=mock_args['nhi_high_cut'])
        script += "&> {path}/merge_dla.log &\n".format(path=mock_args['logs_dir_mergechunks'])
        script += "pid_dla=$!\n"

    if "merge_rand_dla" in todo:
        script += """echo -e "*** Merging DLA randoms ***"\n"""
        if mock_args['use_time']:
            script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
        script += "merge_dla.py -indir {base} -outdir {outpath} -random True ".format(base=mock_args['base_dir'], outpath=mock_args['out_dir'])
        if mock_args['nhi_low_cut'] is not None and mock_args['nhi_high_cut'] is not None:
            script += " --nhi-low-cut {cut1} --nhi-high-cut {cut2} ".format(cut1=mock_args['nhi_low_cut'], cut2=mock_args['nhi_high_cut'])
        script += "&> {path}/merge_rand_dla.log &\n".format(path=mock_args['logs_dir_mergechunks'])
        script += "pid_dla_rand=$!\n"

    if "merge_dla" in todo:
        script += """
if wait $pid_dla; then
    echo "merge_dla OK"
else
    echo "Error in merge_dla"
    exit 1
fi
"""

    if "merge_rand_dla" in todo:
        script += """
if wait $pid_dla_rand; then
    echo "merge_dla randoms OK"
else
    echo "Error in merge_dla randoms"
    exit 1
fi
"""

    if "merge_dla" in todo or "merge_rand_dla" in todo:
        script += """echo -e "==> DLA catalogs done. $(( SECONDS - start )) s"\n"""

    if "transmissions" in todo:
        script += """echo -e "*** Running make_transmissions {threads} times ***"\n""".format(threads=sbatch_args['threads_mergechunks'])
        script += "pids=''\n"
        for job in range(sbatch_args['threads_mergechunks']):
            if mock_args['use_time']:
                script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
            script += "make_transmissions.py -inDir {inpath} -outDir {outpath} -nside {nside} -nest {nest} -job {job} -ncpu {threads} -dla {dla} ".format(inpath=mock_args['base_dir'], outpath=mock_args['out_dir'], nside=mock_args['nside'], nest=mock_args['nest'], job=job, threads=sbatch_args['threads_mergechunks'], dla=mock_args['dla'])
            script += "&> {path}/make_transmissions-{job}.log &\n".format(path=mock_args['logs_dir_mergechunks'], job=job)
            script += """pids+=" $!"\n"""
        script += get_errors("make_transmissions", 0)
        script += """echo -e "==> make_transmissions done. $(( SECONDS - start )) s"\n"""

    script += """echo "----- END -----"\n"""

    filename = mock_args['run_dir']+'/run_mergechunks.sh'
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()


def run_python_script(i_node, i_chunk, codename, mock_args, sbatch_args, name=None):
    """
    write a .sh script that runs codename.py code
    """
    if name is None:
        name = codename
    imin = i_node*sbatch_args['threads_chunk']
    imax = (i_node+1)*sbatch_args['threads_chunk']
    if imax >= mock_args['nslice']: imax=mock_args['nslice']
    script = "#!/bin/bash -l\n"
    script += """pids=""\n"""
    script += """echo "Launching {codename} from {imin} to {imax}"\n""".format(codename=codename, imin=imin, imax=imax-1)
    for job in range(imin, imax):
        if mock_args['use_time']:
            script += """/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " """
        script += "{codename}.py -i {job} ".format(codename=codename, job=job)
        script += mock_args['args_{}'.format(codename)]
        script += " &> {path}/{name}-{job}.log &\n".format(path=mock_args['logs_dir_chunk-'+i_chunk], name=name, job=job)
        script += """pids+=" $!"\n"""
    script += get_errors(codename, imin, threads_num=sbatch_args['threads_chunk'])

    filename = mock_args['run_dir_chunk-'+i_chunk]+'/run_{name}-{i_chunk}-{i_node}.sh'.format(
        name=name, i_chunk=i_chunk, i_node=i_node)
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()


def do_dir_tree(outdir, nside):
    '''
    This function produces the directories to store the transmission files,
    according to the desi format
    '''
    npixel = hp.nside2npix(nside)
    print("Creating directory tree...")
    for i in range(npixel//100 + 1):
        if os.path.isdir(outdir+"/"+str(i)):
            continue
        subprocess.call(["mkdir", outdir+"/"+str(i)])
        # for j in np.unique(pixels[(pixels >= i*100) & (pixels < (i+1)*100)]):
        for j in range(i*100, (i+1)*100):
            if j < npixel:
                subprocess.call(["mkdir", outdir+"/"+str(i)+"/"+str(j)])
    print("Directory tree created.")


def chunk_parameters(cells,stripe_footprint=False):
    '''
    Returns the parameters (ra, dec, ...) of the different chunks according
    to the box size.
    The nominal size is 2560 with 7 chunks. Other size are for debugging.
    '''

    if cells == 2560:
        ra0=['125.5', '189.982649735', '254.46529947', '140', '245.305821271', '-24', '40.4826497351']
        dra=['32.2413248675', '32.2413248675', '32.2413248675', '52.6529106357', '52.6529106357', '32.2413248675', '32.2413248675']
        dec0=['20', '20', '20', '64.4956624338', '64.4956624338', '7', '7']
        ddec=['32.2413248675', '32.2413248675', '32.2413248675', '12.2543375662', '12.2543375662', '32.2413248675', '32.2413248675']
        chunkid=['1', '2', '3', '4', '5', '6', '7']
        nslice = 512

    if cells == 256:
        ra0=['189.982649735', '183.642649735', '177.302649735', '170.962649735', '196.322649735', '202.662649735']
        dra=['3.17', '3.17', '3.17', '3.17', '3.17', '3.17']
        dec0=['20', '20', '20', '20', '20', '20']
        ddec=['3.17', '3.17', '3.17', '3.17', '3.17', '3.17']
        chunkid=['1', '2', '3', '4', '5', '6']
        nslice = 8

    if cells == 512:
        ra0=['190', '202.8', '177.2', '215.6', '190', '202.8', '177.2', '215.6']
        dra=['6.4', '6.4', '6.4', '6.4', '6.4', '6.4', '6.4', '6.4']
        dec0=['0', '0', '0', '0', '12.8', '12.8', '12.8', '12.8']
        ddec=['6.4', '6.4', '6.4', '6.4', '6.4', '6.4', '6.4', '6.4']
        chunkid=['1', '2', '3', '4', '5', '6', '7', '8']
        nslice = 8

    if cells == 1024:
        ra0=['190']
        dra=['12.7']
        dec0=['0']
        ddec=['12.7']
        chunkid=['1']
        nslice = 32

    if cells == 128:
        ra0=['190']
        dra=['1.6']
        dec0=['0']
        ddec=['1.6']
        chunkid=['1']
        nslice = 8

    if cells == 32:
        ra0=['190']
        dra=['1.6']
        dec0=['0']
        ddec=['1.6']
        chunkid=['1']
        nslice = 2

    if stripe_footprint :
        ra0=['-25','25']
        dra=['25','25']
        dec0=['0','0']
        ddec=['2','2']
        chunkid=['1', '2']
        nslice = 512

    return np.array(ra0), np.array(dra), np.array(dec0), np.array(ddec), np.array(chunkid), np.array(nslice)


def make_pk_dir(mock_args):
    '''
    This function creates the directories need for the Pk
    and then writes down the scripts to run the code
    '''
    dir_pk = mock_args['mock_dir']+"/pk"
    try:
        os.makedirs(dir_pk)
    except OSError:
        pass
    mock_args['dir_pk'] = dir_pk
    if mock_args['burst_buffer']:
        mock_args['dir_pk_no_bb'] = dir_pk
    # Create logs directory:
    dir_pk_logs = mock_args['dir_pk']+"/logs"
    try:
        os.makedirs(dir_pk_logs)
    except OSError:
        pass
    mock_args['dir_pk_logs'] = dir_pk_logs
    # Create run directory:
    dir_pk_run = mock_args['dir_pk']+"/run"
    try:
        os.makedirs(dir_pk_run)
    except OSError:
        pass
    mock_args['dir_pk_run'] = dir_pk_run


def make_realisation(imock, mock_args, run_args, sbatch_args):
    '''
    This function creates all the directories needed for a given realisation
    and then writes down the scripts to run the code, as well as the submit.sh
    bash script
    '''
    print("Generating scripts for realisation {}".format(imock))
    mock_args['imock'] = imock
    base_dir = mock_args['mock_dir']+"/mock_{}".format(imock)
    try:
        os.makedirs(base_dir)
    except OSError:
        pass
    mock_args['base_dir'] = base_dir

    ### Define output directories
    if 'special_out' in mock_args.keys():
        if 'out_version' in mock_args.keys():
            out_dir = mock_args['special_out']+'/'+mock_args['out_version']+'.'+imock
        else:
            out_dir = mock_args['special_out']+'/'+imock
        if mock_args['burst_buffer']:
            mock_args['out_dir_no_bb'] = out_dir  # path where the output is copied out of burst buffer
            out_dir = mock_args['base_dir']+"/output"
    else:
        if mock_args['nhi_low_cut'] is not None and mock_args['nhi_high_cut'] is not None:
            out_dir = mock_args['base_dir']+"/output_nhi_{}_{}".format(mock_args['nhi_low_cut'], mock_args['nhi_high_cut'])
        else:
            out_dir = mock_args['base_dir']+"/output"
    try:
        os.makedirs(out_dir)
    except OSError:
        pass
    mock_args['out_dir'] = out_dir

    # desi format directories:
    if run_args['run_mergechunks'] and "transmissions" in run_args['todo_mergechunks']:
        do_dir_tree(out_dir, mock_args['nside'])

    # logs directory:
    logs_dir = out_dir+"/logs"
    try:
        os.makedirs(logs_dir)
    except OSError:
        pass
    mock_args['logs_dir'] = logs_dir

    # runs directory
    run_dir = out_dir+"/runs"
    try:
        os.makedirs(run_dir)
    except OSError:
        pass
    mock_args['run_dir'] = run_dir

    # Chunks directories
    if run_args['run_boxes'] or run_args['run_chunks']:
        for cid in mock_args['chunkid']:
            # Define directories
            chunk_dir = base_dir+"/chunk_"+cid
            dir_boxes = chunk_dir+"/boxes"
            dir_qso = chunk_dir+"/qso/"
            dir_rand = chunk_dir+"/randoms/"
            dir_spectra = chunk_dir+"/spectra"
            dir_spectra_merged = chunk_dir+"/spectra_merged"
            logs_dir_chunk = logs_dir+"/chunk_"+cid
            run_dir_chunk = run_dir+"/chunk_"+cid
            try:
                os.makedirs(chunk_dir)
            except OSError:
                pass
            try:
                os.makedirs(dir_boxes)
            except OSError:
                pass
            try:
                os.makedirs(dir_qso)
            except OSError:
                pass
            try:
                os.makedirs(dir_rand)
            except OSError:
                pass
            try:
                os.makedirs(dir_spectra)
            except OSError:
                pass
            try:
                os.makedirs(dir_spectra_merged)
            except OSError:
                pass
            try:
                os.makedirs(logs_dir_chunk)
            except OSError:
                pass
            try:
                os.makedirs(run_dir_chunk)
            except OSError:
                pass
            mock_args['chunk_dir-'+cid] = chunk_dir
            mock_args['dir_boxes-'+cid] = dir_boxes
            mock_args['dir_qso-'+cid] = dir_qso
            mock_args['dir_rand-'+cid] = dir_rand
            mock_args['dir_spectra-'+cid] = dir_spectra
            mock_args['dir_spectra_merged-'+cid] = dir_spectra_merged
            mock_args['logs_dir_chunk-'+cid] = logs_dir_chunk
            mock_args['run_dir_chunk-'+cid] = run_dir_chunk
    # Mergechunks directories
    if run_args['run_mergechunks']:
        logs_dir_mergechunks = logs_dir+"/mergechunks"
        run_dir_mergechunks = run_dir+"/mergechunks"
        try:
            os.makedirs(logs_dir_mergechunks)
        except OSError:
            pass
        try:
            os.makedirs(run_dir_mergechunks)
        except OSError:
            pass
        mock_args['logs_dir_mergechunks'] = logs_dir_mergechunks
        mock_args['run_dir_mergechunks'] = run_dir_mergechunks

    # Write script for burst buffer
    if mock_args['burst_buffer']:
        i = mock_args['bb_name'].rfind('_')
        if i > 0:
            mock_args['bb_name'] = mock_args['bb_name'][:i]+"_"+imock
        else:
            mock_args['bb_name'] = mock_args['bb_name']+"_"+imock
        # if run_args['run_create']:
        create_reservation(mock_args)
        # if run_args['run_stagein']:
        stage_in(mock_args)
        # if run_args['run_stageout']:
        stage_out(mock_args)
        # if run_args['run_delete']:
        delete_reservation(mock_args)
        # Change the base directory to burst buffer directory:
        for k in mock_args.keys():
            # this dir is not included in the next replace, so we do it by hand
            mock_args['dir_pk'] = "$DW_PERSISTENT_STRIPED_{name}/pk".format(name=mock_args['bb_name'])
            if 'dir' in k:
                if k == 'mock_dir' or 'run' in k or 'log' in k or 'python' in k or 'stage_out' in k or 'no_bb' in k:
                    continue
                mock_args[k] = mock_args[k].replace(mock_args['mock_dir'], "$DW_PERSISTENT_STRIPED_{name}".format(name=mock_args['bb_name']))

    ### Write scripts for each chunks:
    if run_args['run_boxes'] or run_args['run_chunks']:
        for i, cid in enumerate(mock_args['chunkid']):
            mock_args['i_chunk'] = cid
            if run_args['run_boxes']:
                boxes(mock_args, sbatch_args)
            if run_args['run_chunks']:
                chunk(run_args['todo_chunk'], mock_args, sbatch_args)
                for node in range(sbatch_args['nodes_chunk']):
                    if run_args['draw_qso']:
                        mock_args['args_draw_qso'] = "-Nslice "+str(mock_args['nslice'])
                        mock_args['args_draw_qso'] += " -indir "+mock_args['dir_boxes-'+cid]
                        mock_args['args_draw_qso'] += " -outpath "+mock_args['dir_qso-'+cid]
                        mock_args['args_draw_qso'] += " -ra0 "+mock_args['ra0'][i]
                        mock_args['args_draw_qso'] += " -dec0 "+mock_args['dec0'][i]
                        mock_args['args_draw_qso'] += " -chunk "+cid
                        mock_args['args_draw_qso'] += " -dra "+mock_args['dra'][i]
                        mock_args['args_draw_qso'] += " -ddec "+mock_args['ddec'][i]
                        mock_args['args_draw_qso'] += " -zmin "+str(mock_args['zmin'])
                        mock_args['args_draw_qso'] += " -zmax "+str(mock_args['zmax'])
                        mock_args['args_draw_qso'] += " -desi "+str(mock_args['desifootprint'])
                        mock_args['args_draw_qso'] += " -rsd "+str(mock_args['rsd'])
                        # mock_args['args_draw_qso'] += " "+mock_args['seed']+" "+mock_args['zfix']
                        mock_args['args_draw_qso'] += " "+mock_args['seed']
                        run_python_script(node, cid, "draw_qso", mock_args, sbatch_args)
                    if run_args['randoms']:
                        mock_args['args_draw_qso'] = "-Nslice "+str(mock_args['nslice'])
                        mock_args['args_draw_qso'] += " -indir "+mock_args['dir_boxes-'+cid]
                        mock_args['args_draw_qso'] += " -outpath "+mock_args['dir_rand-'+cid]
                        mock_args['args_draw_qso'] += " -ra0 "+mock_args['ra0'][i]
                        mock_args['args_draw_qso'] += " -dec0 "+mock_args['dec0'][i]
                        mock_args['args_draw_qso'] += " -chunk "+cid
                        mock_args['args_draw_qso'] += " -dra "+mock_args['dra'][i]
                        mock_args['args_draw_qso'] += " -ddec "+mock_args['ddec'][i]
                        mock_args['args_draw_qso'] += " -zmin "+str(mock_args['zmin'])
                        mock_args['args_draw_qso'] += " -zmax "+str(mock_args['zmax'])
                        mock_args['args_draw_qso'] += " -desi "+str(mock_args['desifootprint'])
                        mock_args['args_draw_qso'] += " -rsd "+str(mock_args['rsd'])
                        # mock_args['args_draw_qso'] += " "+mock_args['seed']+" "+mock_args['zfix']
                        mock_args['args_draw_qso'] += " "+mock_args['seed']
                        mock_args['args_draw_qso'] += " -random True "
                        run_python_script(node, cid, "draw_qso", mock_args, sbatch_args, "randoms")
                    if run_args['make_spectra']:
                        mock_args['args_make_spectra'] = "-QSOfile "+mock_args['dir_qso-'+cid]+"/QSO-"
                        mock_args['args_make_spectra'] += " -boxdir "+mock_args['dir_boxes-'+cid]
                        mock_args['args_make_spectra'] += " -outDir "+mock_args['dir_spectra-'+cid]
                        mock_args['args_make_spectra'] += " -N "+str(mock_args['nslice'])
                        mock_args['args_make_spectra'] += " -zmin "+str(mock_args['zmin'])
                        mock_args['args_make_spectra'] += " -zmax "+str(mock_args['zmax'])
                        mock_args['args_make_spectra'] += " -rsd "+str(mock_args['rsd'])
                        mock_args['args_make_spectra'] += " -dla "+str(mock_args['dla'])
                        run_python_script(node, cid, "make_spectra", mock_args, sbatch_args)
                    if run_args['merge_spectra']:
                        mock_args['args_merge_spectra'] = "-inDir "+mock_args['dir_spectra-'+cid]
                        mock_args['args_merge_spectra'] += " -outDir "+mock_args['dir_spectra_merged-'+cid]
                        mock_args['args_merge_spectra'] += " -aa "+str(mock_args['a'])
                        mock_args['args_merge_spectra'] += " -bb "+str(mock_args['b'])
                        mock_args['args_merge_spectra'] += " -cc "+str(mock_args['c'])
                        mock_args['args_merge_spectra'] += " -nside "+str(mock_args['nside'])
                        mock_args['args_merge_spectra'] += " -nest "+str(mock_args['nest'])
                        mock_args['args_merge_spectra'] += " -rsd "+str(mock_args['rsd'])
                        mock_args['args_merge_spectra'] += " -addnoise "+str(mock_args['small_scales'])
                        mock_args['args_merge_spectra'] += " -dla "+str(mock_args['dla'])
                        mock_args['args_merge_spectra'] += " --store-g "+str(mock_args['store_g'])
                        if mock_args['p1dfile'] is not None:
                            mock_args['args_merge_spectra'] += " --fit-p1d True "  # in order to read the p1dfile format
                            mock_args['args_merge_spectra'] += " -p1dfile "+mock_args['p1dfile']                            
                        if mock_args['fit_p1d']:
                            mock_args['args_merge_spectra'] += " --fit-p1d True "
                            mock_args['args_merge_spectra'] += " -p1dfile "+mock_args['chunk_dir-'+cid]+"/p1dmiss.fits"
                        mock_args['args_merge_spectra'] += " "+mock_args['seed']+" "+mock_args['zfix']
                        run_python_script(node, cid, "merge_spectra", mock_args, sbatch_args)
    if run_args['run_mergechunks']:
        mergechunks(run_args['todo_mergechunks'], mock_args, sbatch_args)
    submit(mock_args, run_args)
    print("{}/submit.sh generated. Run it to submit the sbatch jobs of this realisation to the cori queue\n".format(mock_args['run_dir']))
    # The pk needs to be run only once:
    run_args['run_pk'] = False


def submit(mock_args, run_args):
    '''
    This function write the bash script submit.sh, which is used to submit
    the different sbatch script to cori nodes
    '''
    path = mock_args['run_dir']
    script = "#!/bin/bash\n"
    if mock_args['sbatch']:
        # Pk
        if run_args['run_pk']:
            script += "run_pk=$(sbatch --parsable "
            script += "--output "+mock_args['dir_pk_logs']+"/run_pk.log "
            script += mock_args['dir_pk_run']+"/run_pk.sh)\n"
            script += """echo "run_pk.sh: "$run_pk\n"""
        # run create
        if mock_args['burst_buffer'] and run_args['run_create']:
            script += "run_create=$(sbatch --parsable "
            if run_args['run_pk']:
                script += "--dependency=afterok:$run_pk "
            script += "--output "+mock_args['logs_dir']+"/run_create.log "
            script += path+"/create_reservation.sh)\n"
            script += """echo "run_create.sh: "$run_create \n"""
        # run stage in
        if mock_args['burst_buffer'] and run_args['run_stagein']:
            script += "run_stagein=$(sbatch --parsable "
            script += "--output "+mock_args['logs_dir']+"/run_stagein.log "
            if  run_args['run_create'] and mock_args['burst_buffer']:
                script += "--dependency=afterok:$run_create "
            elif run_args['run_pk']:
                script += "--dependency=afterok:$run_pk "
            script += path+"/stage_in.sh)\n"
            script += """echo "stage_in.sh: "$run_stagein \n"""
        # run boxes and chunks
        for i, cid in enumerate(mock_args['chunkid']):
            if run_args['run_boxes']:
                script += "run_boxes_{i}=$(sbatch --parsable ".format(i=cid)
                if run_args['run_pk'] or ((run_args['run_create'] or run_args['run_stagein']) and mock_args['burst_buffer']):
                    script += "-d afterok:"
                    afterok = ""
                    if run_args['run_pk']:
                        afterok += "$run_pk,"
                    if mock_args['burst_buffer']:
                        if run_args['run_create']:
                            afterok += "$run_create,"
                        if run_args['run_stagein']:
                            afterok += '$run_stagein,'
                    script += afterok[:-1]
                script += " --output "+mock_args['logs_dir']+"/run_boxes-{i}.log ".format(i=cid)
                script += path+"/run_boxes-{i}.sh)\n".format(i=cid)
                script += """echo "run_boxes-{i}.sh: "$run_boxes_{i}\n""".format(i=cid)
            if run_args['run_chunks']:
                script += "run_chunk_{i}=$(sbatch --parsable ".format(i=cid)
                if run_args['run_boxes'] or ((run_args['run_create'] or run_args['run_stagein']) and mock_args['burst_buffer']):
                    script += "-d afterok:"
                    afterok = ""
                    if run_args['run_boxes']:
                        afterok += "$run_boxes_{i},".format(i=cid)
                    if mock_args['burst_buffer']:
                        if run_args['run_create']:
                            afterok += "$run_create,"
                        if run_args['run_stagein']:
                            afterok += "$run_stagein,"
                    script += afterok[:-1]
                script += " --output "+mock_args['logs_dir']+"/run_chunks-{i}.log ".format(i=cid)
                script += path+"/run_chunk-{i}.sh)\n".format(i=cid)
                script += """echo "run_chunk-{i}.sh: "$run_chunk_{i}\n""".format(i=cid)
        # run mergechunks
        if run_args['run_mergechunks']:
            script += "run_mergechunks=$(sbatch --parsable "
            if run_args['run_chunks'] or ((run_args['run_create'] or run_args['run_stagein']) and mock_args['burst_buffer']):
                script += "-d afterok:"
                afterok = ""
                if run_args['run_chunks']:
                    for cid in mock_args['chunkid']:
                        afterok += "$run_chunk_{i},".format(i=cid)
                if mock_args['burst_buffer']:
                    if run_args['run_create']:
                        afterok += "$run_create,"
                    if run_args['run_stagein']:
                        afterok += "$run_stagein,"
                script += afterok[:-1]
            script += " --output "+mock_args['logs_dir']+"/run_mergechunks.log "
            script += path+"/run_mergechunks.sh)\n"
            script += """echo "run_mergechunks.sh: "$run_mergechunks\n"""
        # Burst buffer jobs:
        if mock_args['burst_buffer'] and run_args['run_stageout']:
            script += "run_stageout=$(sbatch --parsable "
            if run_args['run_boxes'] or run_args['run_chunks'] or run_args['run_mergechunks']:
                script += "-d afterany:"
                afterany = ""
                for cid in mock_args['chunkid']:
                    if run_args['run_boxes']:
                        afterany += "$run_boxes_{i},".format(i=cid)
                    if run_args['run_chunks']:
                        afterany += "$run_chunk_{i},".format(i=cid)
                if run_args['run_mergechunks']:
                    afterany += "$run_mergechunks "
                script += afterany[:-1]
            script += " --output "+mock_args['logs_dir']+"/run_stageout.log "
            script += path+"/stage_out.sh)\n"
            script += """echo "run_stageout.sh: "$run_stageout \n"""
        if mock_args['burst_buffer'] and run_args['run_delete']:
            script += "run_delete=$(sbatch --parsable "
            if run_args['run_boxes'] or run_args['run_chunks'] or run_args['run_mergechunks'] or run_args['run_stageout']:
                script += "-d afterok:"
                afterok = ""
                for cid in mock_args['chunkid']:
                    if run_args['run_boxes']:
                        afterok += "$run_boxes_{i},".format(i=cid)
                    if run_args['run_chunks']:
                        afterok += "$run_chunk_{i},".format(i=cid)
                if run_args['run_mergechunks']:
                    afterok += "$run_mergechunks,"
                if run_args['run_stageout']:
                    afterok += "$run_stageout "
                script += afterok[:-1]
            script += " --output "+mock_args['logs_dir']+"/run_delete.log "
            script += path+"/delete_reservation.sh)\n"
            script += """echo "run_delete.sh: "$run_delete \n"""
    else:
        if run_args['run_pk']:
            script += "bash {}".format(mock_args['dir_pk_run']+"/run_pk.sh")
            script += " 2>&1 | tee "+mock_args['dir_pk_logs']+"/run_pk.log\n"
        for cid in mock_args['chunkid']:
            if run_args['run_boxes']:
                script += "bash {path}/run_boxes-{i}.sh".format(path=path, i=cid)
                script += " 2>&1 | tee "+mock_args['logs_dir']+"/run_boxes-{i}.log\n".format(i=cid)
            if run_args['run_chunks']:
                script += "bash {path}/run_chunk-{i}.sh".format(path=path, i=cid)
                script += " 2>&1 | tee "+mock_args['logs_dir']+"/run_chunks-{i}.log\n".format(i=cid)
        if run_args['run_mergechunks']:
            script += "bash {}".format(path+"/run_mergechunks.sh")
            script += " 2>&1 | tee "+mock_args['logs_dir']+"/run_mergechunks.log\n"
    filename = mock_args['run_dir']+'/submit.sh'
    fout = open(filename, 'w')
    fout.write(script)
    fout.close()
    os.chmod(filename,os.path.stat.S_IRWXU | os.path.stat.S_IRWXG)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Writes scripts to send Saclay mocks production.\n'+\
                        'It is designe to work on nersc cori nodes.\n'+\
                        'Once done, run submit.sh')

    parser.add_argument("--mock-dir", type=str, default=None, required=True,
        help="directory where all files are saved")

    parser.add_argument("--cori-nodes", type=str, default='True', required=False,
        help="If True, send the code to cori nodes")

    parser.add_argument("--out-dir", type=str, default=None, required=False,
        help="directory where the output in desi format is saved (optional). Default is within the base directory")

    parser.add_argument("--out-version", type=str, default=None, required=False,
        help="Prefix of output directories for the different realisations (optional)")

    parser.add_argument("--box-size", type=int, default=2560, required=False,
        help="Number of cell that defines the chunk box (optional)")

    parser.add_argument("--chunk-id", type=int, nargs="*", default=None, required=False,
        help="Specify the chunks to be produced (optional)")

    parser.add_argument("--realisation-number", type=int, default=1, required=False,
        help="The number of realisation to be produced (optional)")

    parser.add_argument("--realisation-id", nargs="*", default=None, required=False,
        help="Specify particular realisations to be produced (optional)")

    parser.add_argument("--stage-out", type=str, nargs="*", default='all', required=False,
        help="Specify what output you want to save to your scratch directory when using the burst buffer (optional)\n"+
                        "you can backup among: 'boxes qso randoms spectra spectra_merged dla dla_randoms' or just 'all'")

    parser.add_argument("--discard-boxes", action='store_true', required=False,
        help="Do not stage_out the boxes (stage_out only boxk.npy); only usable in BurstBuffer mode")

    parser.add_argument("--bb-nodes", action='store_true', required=False,
        help="use the BurstBuffer nodes")

    parser.add_argument("--fit-p1d", type=str, nargs="*", default = '', required=False,
        help="Set up parameters for the 1D power spectrum fitting procedure:\n"
            +"give in order: z, a, b, c, seed")

    parser.add_argument("--p1d-file", type=str, default=None, required=False,
        help="p1dfile to generate small scales. Default is to read from etc/pkmiss_interp.fits.gz")

    parser.add_argument("--seed", type=int, default=None, required=False,
        help="Specify a particular seed (optional)")

    parser.add_argument("--time", type=str, default="True", required=False,
        help="If True, use /usr/bin/time to time the jobs (optional)")

    parser.add_argument("--account", type=str, default="desi", required=False,
        help="account to be used on cori nodes (optional)")

    parser.add_argument("--email", type=str, default=None, required=False,
        help="Your email address (optional)")

    parser.add_argument("--stripe82", action='store_true', required=False,
        help="Create the box limits adapted to the Stripe 82 footprint")

    args = parser.parse_args()

    ### Define sbatch parameters :
    # These parameters are set up for nominal chunk size (2560)
    sbatch_args = {}
    sbatch_args['account'] = args.account
    sbatch_args['email'] = args.email
    # Parameters for pk job:
    sbatch_args['time_pk'] = "00:20:00"  # default "00:15:00"
    sbatch_args['queue_pk'] = "debug"  # default "regular"
    sbatch_args['name_pk'] = "pk_saclay"
    sbatch_args['threads_pk'] = 16  # default 16
    sbatch_args['nodes_pk'] = 1  # default 1
    # Parameters for box jobs:
    sbatch_args['time_boxes'] = "02:00:00"  # default "01:30:00"
    sbatch_args['queue_boxes'] = "regular"  # default "regular"
    sbatch_args['name_boxes'] = "boxes_saclay"
    sbatch_args['threads_boxes'] = 64  # default 64
    sbatch_args['nodes_boxes'] = 1  # default 1
    # Parameters for chunk jobs:
    sbatch_args['time_chunk'] = "00:30:00"  # default "00:30:00"
    sbatch_args['queue_chunk'] = "regular"  # default "regular"
    sbatch_args['name_chunk'] = "chunk_saclay"
    sbatch_args['threads_chunk'] = 32  # default 32
    sbatch_args['nodes_chunk'] = 16  # nodes * threads should be = nslice, default 16
    # Parameters for mergechunks job:
    sbatch_args['time_mergechunks'] = "01:30:00"  # default "01:30:00"
    sbatch_args['queue_mergechunks'] = "regular"  # default "regular"
    sbatch_args['name_mergechunks'] = "mergechunks_saclay"
    sbatch_args['threads_mergechunks'] = 64  # default 64
    sbatch_args['nodes_mergechunks'] = 1  # default 1

    ### Define mock parameters:
    mock_args = {}
    # Healpix:
    mock_args['nside'] = 16  # healpix nside parameter, default 16
    mock_args['nest'] = True  # healpix nest parameter, default True
    # Box and spectrum parameters:
    mock_args['nx'] = args.box_size
    mock_args['ny'] = args.box_size
    mock_args['nz'] = 1536  # This one is fixed to get the whole redshift range
    mock_args['pixel_size'] = 2.19  # pixel size in Mpc/h
    mock_args['a'] = -1  # a parameter in FGPA; -1 is to read a(z) from etc/params.fits
    mock_args['b'] = -1  # b parameter in FGPA; -1 is to read b(z) from etc/params.fits
    mock_args['c'] = -1  # c paramter in FGPA; -1 is to read c(z) from etc/params.fits
    mock_args['zmin'] = 1.8  # minimal redshift to draw QSO
    mock_args['zmax'] = 3.6  # maximal redshift to draw QSO
    mock_args['zfix'] = "-zfix 2.4"  # "-zfix 2.6" to fix the redshift to a special value; and "" for none
    mock_args['store_g'] = False  # If true, store delta_l, delta_s and eta_par fields
    # mock options:
    mock_args['seed'] = ""  # "-seed 10" to specify a seed, "" to specify nothing
    mock_args['desifootprint'] = True  # If True, cut QSO outside desi footprint
    mock_args['NQSO'] = -1  # If >0, limit the number of QSO treated in make_spectra
    mock_args['small_scales'] = True  # If True, add small scales in FGPA
    mock_args['rsd'] = True  # If True, add RSD
    mock_args['dla'] = True  # If True, add DLA
    mock_args['nmin'] = 17.2  # log(N_HI) min for DLA
    mock_args['nmax'] = 22.5  # log(N_HI) max for DLA
    mock_args['nhi_low_cut'] = None  # cut DLAs with log(n_HI) < cut
    mock_args['nhi_high_cut'] = None  # cut DLAs with log(n_HI) > cut
    # Run options:
    mock_args['use_time'] = util.str2bool(args.time)  # If True, use /usr/bin/time/ to time jobs
    mock_args['verbosity'] = None  # Set it to "-v -v -v -v" if you want info from sbatch jobs
    mock_args['sbatch'] = util.str2bool(args.cori_nodes)  # If True, jobs are sent to cori nodes (frontend nodes otherwise)
    # Burst buffer options:
    mock_args['burst_buffer'] = args.bb_nodes  # If True, use the burst buffer on cori nodes. /!\ only if sbatch is True
    mock_args['bb_size'] = "5TB"  # A mock realisation at nominal size is 4Tb, so ask for 5
    mock_args['bb_name'] = "saclaymock"  # Each realisation has a reservation named 'bb_name-'+i_realisation

    ### Code to runs:
    run_args = {}
    # pk:
    run_args['run_pk'] = False  # Produce Pk
    # boxes:
    run_args['run_boxes'] = False  # Produce GRF boxes
    # chunks:
    run_args['run_chunks'] = False  # produce chunks
    run_args['draw_qso'] = False  # run draw_qso.py
    run_args['randoms'] = False  # run draw_qso.py for randoms
    run_args['make_spectra'] = False  # run make_spectra.py
    run_args['merge_spectra'] = True  # run merge_spectra.py
    # merge chunks:
    run_args['run_mergechunks'] = True  # Gather outputs from all chunks and write in desi format
    run_args['merge_qso'] = True  # Compute master.fits file
    run_args['merge_randoms'] = False  # Compute master_randoms.fits file
    run_args['compute_dla'] = True  # Compute dla catalog of each chunks
    run_args['dla_randoms'] = False  # Compute dla randoms catalogs of each chunks
    run_args['merge_dla'] = True  # Compute master_DLA.fits file
    run_args['merge_rand_dla'] = False  # Compute master_DLA_randoms.fits file
    run_args['transmissions'] = True  # Write transmissions files
    # burst buffer
    run_args['run_create'] = True  # Create persistent reservation
    run_args['run_stagein'] = True  # Stage in the init files (pk, directories, ...) (from scratch to BB)
    run_args['run_stageout'] = True  # Stage out the produced files (from BB to scratch)
    run_args['run_delete'] = False  # delete the persistent reservation

    # -------------------------- Nothing to change bellow
    ### Define directories
    nmocks = args.realisation_number
    if nmocks > 1:
        print("WARNING: Make sure that the Pk run is completed before sending"\
              +" the run of next realisations.\n The Pk will be produce with"\
              +" the run of the first realisation.")
    mock_dir = os.path.expandvars(args.mock_dir)
    mock_args['nmocks'] = nmocks
    mock_args['mock_dir'] = mock_dir

    ### Define chunks parameters :
    if args.seed is not None:
        mock_args['seed'] = "-seed "+str(args.seed)
    if args.box_size < 2560:
        mock_args['desifootprint'] = False
    ra0, dra, dec0, ddec, chunkid, nslice = chunk_parameters(args.box_size,stripe_footprint=args.stripe82)
    if args.chunk_id is not None:
        m = np.array(args.chunk_id) - 1
        ra0 = ra0[m]
        dec0 = dec0[m]
        dra = dra[m]
        ddec = ddec[m]
        chunkid = chunkid[m]

    ### 1D power spectrum fitting procedure
    mock_args['p1dfile'] = args.p1d_file
    mock_args['fit_p1d'] = False
    if args.fit_p1d:
        if len(args.fit_p1d) != 5:
            raise Exception("Number of parameters given with --fit-p1d not correct")
        print("P1D fitting procedure")
        print("zfix, a, b, c, seed = {}".format(args.fit_p1d))
        mock_args['fit_p1d'] = True
        mock_args['nx'] = 256
        mock_args['ny'] = 256
        mock_args['nz'] = 1536
        mock_args['zfix'] = "-zfix {}".format(args.fit_p1d[0])
        mock_args['a'] = args.fit_p1d[1]
        mock_args['b'] = args.fit_p1d[2]
        mock_args['c'] = args.fit_p1d[3]
        mock_args['seed'] = '-seed {}'.format(args.fit_p1d[4])
        mock_args['desifootprint'] = False
        mock_args['dla'] = False
        mock_args['use_time'] = False
        ra0, dra, dec0, ddec, chunkid, nslice = chunk_parameters(mock_args['nx'])
        ra0 = ra0[[0]]
        dra = dra[[0]]
        dec0 = dec0[[0]]
        ddec = ddec[[0]]
        chunkid = chunkid[[0]]
        mock_args['sbatch'] = False
        run_args['run_pk'] = True
        run_args['run_boxes'] = True
        run_args['run_chunks'] = True
        run_args['draw_qso'] = True
        run_args['randoms'] = False
        run_args['make_spectra'] = True
        run_args['merge_spectra'] = True
        run_args['run_mergechunks'] = False

    ### Define sbatch options
    if not mock_args['sbatch']:
        print("Warning: the jobs will not be sent to cori nodes, they will be executed here.")
        mock_args['burst_buffer'] = False
        for k in sbatch_args.keys():
            if 'nodes' in k:
                sbatch_args[k] = 1

    if mock_args['burst_buffer']:
        print("WARNING: You are using Burst Buffer mode. If you cancel the sbatch"\
              +" jobs or if there is any job bug, please cancel the burst buffer"\
              +" persistent reservation by running the script delete_reservation.sh")
        if 'cscratch1' not in args.mock_dir:
            print("Error: --mock-dir option should point to /global/cscratch1 when using the burst buffer mode !")
            sys.exit(1)
        mock_args['stage_out_dir'] = args.stage_out
        mock_args['discard_boxes'] = args.discard_boxes
        if mock_args['discard_boxes']:
            print("WARNING: boxes will be discarded.")
        if ' ' in mock_args['stage_out_dir'][0] or ',' in mock_args['stage_out_dir'][0]:
            print("Error: --stage-out option should be given as: --stage-out boxes qso ... (without global quotes)")
            sys.exit(1)
    else:
        run_args['run_create'] = False
        run_args['run_stageout'] = False
        run_args['run_delete'] = False

    mock_args['ra0'] = ra0
    mock_args['dec0'] = dec0
    mock_args['dra'] = dra
    mock_args['ddec'] = ddec
    mock_args['chunkid'] = chunkid
    mock_args['nslice'] = nslice

    ### Codes to run
    run_args['todo_chunk'] = ""
    if run_args['draw_qso']: run_args['todo_chunk'] += "qso "
    if run_args['randoms']: run_args['todo_chunk'] += "randoms "
    if run_args['make_spectra']: run_args['todo_chunk'] += "make_spectra "
    if run_args['merge_spectra']: run_args['todo_chunk'] += "merge_spectra "

    run_args['todo_mergechunks'] = ""
    if run_args['merge_qso']: run_args['todo_mergechunks'] += "merge_qso "
    if run_args['merge_randoms']: run_args['todo_mergechunks'] += "merge_randoms "
    if run_args['compute_dla'] and mock_args['dla']: run_args['todo_mergechunks'] += "compute_dla "
    if run_args['dla_randoms'] and mock_args['dla']: run_args['todo_mergechunks'] += "dla_randoms "
    if run_args['merge_dla'] and mock_args['dla']: run_args['todo_mergechunks'] += "merge_dla "
    if run_args['merge_rand_dla'] and mock_args['dla']: run_args['todo_mergechunks'] += "merge_rand_dla "
    if run_args['transmissions']: run_args['todo_mergechunks'] += "transmissions "

    print("Writting scripts for {} realisations".format(nmocks))
    print("Mock files will be written in {}".format(mock_dir))
    if args.out_dir is not None:
        mock_args['special_out'] = os.path.expandvars(args.out_dir)
        if args.out_version is not None:
            mock_args['out_version'] = os.path.expandvars(args.out_version)
        print("Outputs will be written in {}".format(mock_args['special_out']))
    print("It will run:")
    if run_args['run_pk']: print("- Pk")
    if mock_args['burst_buffer'] and run_args['run_create']: print("- run_create")
    if mock_args['burst_buffer'] and run_args['run_stagein']: print("- run_stagein")
    if run_args['run_boxes']: print("- Boxes")
    if run_args['run_chunks']: print("- Chunks: {}".format(run_args['todo_chunk']))
    if run_args['run_mergechunks']: print("- Merge Chunks: {}".format(run_args['todo_mergechunks']))
    if mock_args['burst_buffer'] and run_args['run_stageout']: print("- Stage out: {}".format(mock_args['stage_out_dir']))
    if mock_args['burst_buffer'] and run_args['run_delete']: print("- Delete reservation")

    print("It will produce {} chunks, with ids: {}".format(len(chunkid), chunkid))
    print("and ra0: {}".format(ra0))

    # define and create directories for Pk:
    make_pk_dir(mock_args)
    # Write Pk scripts:
    if run_args['run_pk']:
        print("Generating script for Pk...\n")
        pk(mock_args, sbatch_args)

    # Write scripts for all realisations
    if args.realisation_id is not None:
        for imock in args.realisation_id:
            make_realisation(imock, mock_args, run_args, sbatch_args)
    else:
        for imock in range(nmocks):  # loop over realisations
            make_realisation(imock, mock_args, run_args, sbatch_args)


if __name__ == "__main__":
    main()
