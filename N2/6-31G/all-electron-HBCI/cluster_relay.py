"""Run on the workstation with authenticated SSH access to both Sparks."""
import argparse
import fcntl
import json
from pathlib import Path
import subprocess
import time

PRIMARY, SECONDARY = "spark-66ea", "spark-6ad0"
BASE = "MBRefData-runs/N2/6-31G/"
QUEUE = BASE + "all-electron-HBCI/cluster_queue.py"
SSH = ["ssh","-o","BatchMode=yes","-o","ConnectTimeout=10","-o","ServerAliveInterval=10",
       "-o","ServerAliveCountMax=2"]


def call(host, *args):
    return json.loads(subprocess.check_output([*SSH,host,"python3",QUEUE,*args], timeout=60))


def transfer(source, dest, paths, checksum=False):
    if paths:
        cmd = ["rsync","-a","--files-from=-","-e"," ".join(SSH)]
        if checksum:
            cmd.append("--checksum")
        subprocess.run([*cmd,source,dest], input=("\n".join(paths)+"\n").encode(),
                       timeout=240, check=True)


def cycle(cache):
    remote1, remote2 = PRIMARY+":"+BASE, SECONDARY+":"+BASE
    local = str(cache)+"/"
    second = call(SECONDARY,"tick")
    while second["free_slots"]:
        r = call(PRIMARY,"reserve",*second["owners"])
        if not r:
            break
        seed = call(PRIMARY,"seed",r)
        transfer(remote1,local,seed["paths"],True)
        transfer(local,remote2,seed["paths"],True)
        print("Verified seed",call(SECONDARY,"ready",r,seed["sha"]),flush=True)
        second = call(SECONDARY,"tick")
        print("Dispatched",second["launched"],flush=True)
    exports = call(SECONDARY,"export")
    # Fit files change in place. Immutable raw results and logs go first.
    for key in ("raw","fits"):
        transfer(remote2,local,exports[key],key=="fits")
        transfer(local,remote1,exports[key],key=="fits")
    first = call(PRIMARY,"status")
    if first["failures"] or second["failures"]:
        print("Failed sectors", first["failures"],second["failures"],flush=True)
    return first["complete"]


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    args=parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    with (args.cache/"broker.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        while True:
            try:
                if cycle(args.cache.resolve()):
                    break
            except (OSError, ValueError, subprocess.SubprocessError) as e:
                print(time.ctime(),e,flush=True)
            time.sleep(45)
