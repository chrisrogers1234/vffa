import time
import os
import sys
import subprocess
import glob

PROC_QUEUE = []
PROC_RUNNING = []
UNIQUE_ID = 0
N_PROCS = 5
TARGET_SCRIPT = "run_sim.py"
TIME_0 = time.time()

def will_make_new_procs(temp_proc_queue):
    global PROC_QUEUE, N_PROCS
    return len(temp_proc_queue) < N_PROCS and len(PROC_QUEUE) > 0

def poll_process_queue():
    global UNIQUE_ID, PROC_RUNNING, PROC_QUEUE, TIME
    if is_scarf():
        temp_proc_queue = poll_scarf()
    else:
        temp_proc_queue = poll_laptop()
    print("\r", round(time.time()-TIME_0, 1), "...",
          "Running", len(PROC_RUNNING), 
          "with", len(PROC_QUEUE), "queued", end=" ")
    if will_make_new_procs(temp_proc_queue):
        print()
    while will_make_new_procs(temp_proc_queue):
        subproc_args, logname = PROC_QUEUE.pop(0)
        UNIQUE_ID += 1
        job_log = "logs/"+logname+".log"
        if is_scarf():
            subproc_args = ['bsub',
                            '-q', 'scarf-ibis',
                            '-n', '1',
                            '-W', '72:00',
                            '-o', job_log]+subproc_args
            logfile = open("logs/"+logname+".bsub", "w")
        else:
            logfile = open(job_log, "w")
        temp_proc_queue.append(subprocess.Popen(subproc_args,
                               stdout=logfile, stderr=subprocess.STDOUT))
        print("Running", subproc_args, "with log", job_log,
              "pid", temp_proc_queue[-1].pid, len(PROC_RUNNING))
    PROC_RUNNING = temp_proc_queue

def poll_laptop():
    temp_proc_queue = []
    for proc in PROC_RUNNING:
        if proc.poll() == None:
            temp_proc_queue.append(proc)
        else:
            print("\nPID", proc.pid, "finished with return code", proc.returncode)
    sys.stdout.flush()
    return temp_proc_queue

def poll_scarf():
    global TARGET_SCRIPT
    output = subprocess.check_output(['bjobs', '-w'])
    lines = output.split('\n')[1:]
    script_lines = []
    for line in lines:
        if TARGET_SCRIPT not in line:
            continue
        p_index = line.index('python')
        line = line[p_index:]
        line.split(' ')
        script_lines.append(line)
    return script_lines

def is_scarf():
    uname = str(subprocess.check_output(['uname', '-a']))
    return uname.find('scarf.rl.ac.uk') > -1

def main(configs):
    if os.getenv("OPAL_EXE_PATH") == None:
        raise ValueError("No OPAL_EXE_PATH set")
    global N_PROCS, TARGET_SCRIPT, UNIQUE_ID
    if not os.path.exists("logs"):
        os.makedirs("logs")
    if is_scarf():
        N_PROCS = 150
    for config in configs:
        log_file = config.split("/")[-1]
        log_file = log_file[:-3]
        proc_tuple = (["python", "scripts/bin/run_one.py", config], log_file)
        PROC_QUEUE.append(proc_tuple)
    print(len(PROC_QUEUE), "jobs")
    while len(PROC_QUEUE) > 0 or len(PROC_RUNNING) > 0:
        poll_process_queue()
        if len(PROC_QUEUE) == 0 and len(PROC_RUNNING) == 0:
            break
        time.sleep(5)

if __name__ == "__main__":
    main(sys.argv[1:])
    print("\nFinished")



