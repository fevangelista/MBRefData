#!/usr/bin/env python3
"""Read-only terminal monitor for the N₂ HBCI campaign on both Spark hosts."""
import argparse
import concurrent.futures
import curses
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HOSTS = ("spark-66ea", "spark-6ad0")
SECTORS = ("singlet-Ag", "triplet-B1u")
LABELS = ("¹Ag", "³B₁ᵤ")
DATASET = Path("N2/6-31G/all-electron-HBCI")
SCRIPT = Path(__file__).resolve()


def read_json(path, default=None):
    try:
        return json.loads(path.read_text())
    except FileNotFoundError:
        return default


def host_grid(root, host):
    grid = read_json(root.parent / "all-electron-FCI/input.json")
    return grid["distances_angstrom" if host == HOSTS[0] else "additional_hbci_distances_angstrom"]


def tail(path, size=8000):
    try:
        with path.open("rb") as stream:
            stream.seek(max(0, path.stat().st_size - size))
            return stream.read().decode(errors="replace").splitlines()
    except FileNotFoundError:
        return []


def solver_stage(path):
    for line in reversed(tail(path)):
        row = re.match(r"\s*(\d+)\s+(-\d+\.\d+)\s+[-+\d.eE]+\s+([-+\d.eE]+)\s+\d+\s+[\d.]+", line)
        if row:
            return f"Davidson iter {row[1]}   residual {row[3]} Eh"
        if "processed" in line and "batches" in line:
            return "Scanning external determinants"
        if "macroiteration" in line.lower() or "macro iteration" in line.lower():
            return line.strip()[:90]
        if line.startswith("Finished "):
            return "Writing completed threshold"
    return "Starting calculation"


def jobs_on_host(root):
    jobs = []
    clock = os.sysconf("SC_CLK_TCK")
    uptime = float(Path("/proc/uptime").read_text().split()[0])
    for proc in Path("/proc").glob("[0-9]*"):
        try:
            args = (proc/"cmdline").read_bytes().decode(errors="replace").split("\0")
            scripts = {"threshold.py", "precision_fci.py", "singlet_spectrum.py", "calibrate_fci.py"}
            i = next((i for i,a in enumerate(args) if a.rsplit("/",1)[-1] in scripts), None)
            if i is None or not Path(args[0]).name.startswith("python"):
                continue
            script=Path(args[i])
            if not script.is_absolute():
                script=Path(os.readlink(proc/"cwd"))/script
            if script.parent not in (root,root/"validation"):
                continue
            eps=None
            buffered=False
            if script.name=="threshold.py":
                space, r, sector, raw_eps = args[i+1:i+5]
                eps=float(raw_eps)
                buffered = "--buffer-roots" in args
                folder = root/space/f"R-{float(r):.2f}"
                if buffered:
                    folder = root/"diagnostics"/folder.name if space == "full" else folder/"buffered"
                name = f"{sector}-{eps:.8g}" + ("-pt2" if "--tight-pt2-copy" in args else "")
                log = folder/(name+".out")
                kind=("full" if space=="full" else "reduced")+" "+("4 roots" if buffered else "2 roots")
            elif script.name=="singlet_spectrum.py":
                space,r="validation",4.0
                sector="triplet-B1u" if "--triplet" in args else "singlet-Ag"
                log=root/space/"R-4.00"/("triplet-spectrum.out" if "--triplet" in args else "singlet-spectrum.out")
                kind="exact FCI 4 roots"
            else:
                if script.name=="precision_fci.py":
                    space,r=args[i+1:i+3]
                    filename="fci-tight.out"
                else:
                    r=args[i+1]
                    space="validation" if "--validation" in args else "calibration"
                    filename="fci.out"
                sector="both sectors"
                log=root/space/f"R-{float(r):.2f}"/filename
                kind="exact FCI control"
            # Another coordinator can wait on the same calculation's file lock.
            if "lock" in (proc/"wchan").read_text():
                continue
            stat = (proc/"stat").read_text().rsplit(")", 1)[1].split()
            jobs.append(dict(pid=int(proc.name), r=float(r), sector=sector, space=space,
                             buffer=buffered, eps=eps, kind=kind, stage=solver_stage(log),
                             elapsed=max(0,uptime-int(stat[19])/clock), log=str(log)))
        except (OSError, ValueError, IndexError):
            continue  # A process can exit between reads.
    return jobs


def rank_checks(root, distance):
    point = root/"diagnostics"/f"R-{distance:.2f}"
    checks = {}
    for name in ("initial-comparison.json", "results.json", "triplet-rank-check.json"):
        check = read_json(point/name, {})
        if "rank_confirmed" in check:
            checks[check["sector"]] = check
    return checks


def sector_status(root, distance, sector, target, jobs, check, failures):
    point = root/"full"/f"R-{distance:.2f}"
    fit = read_json(point/(sector+"-fit.json"), {})
    current = fit.get("target_Eh") == target
    fitted = fit.get("fits", []) if current else []
    uncertainty = [s["estimated_uncertainty_Eh"]*1e6 for s in fitted] or [None, None]
    passed = [bool(fit.get("qualified") and s["diagnostics_pass"] and s["estimated_uncertainty_Eh"] <= target)
              for s in fitted] or [False, False]
    verified_rank = None
    if check:
        verified_rank = bool(check.get("rank_confirmed") and check.get("target_Eh") == target)
        j = check["rank"]-1
        passed[j] = verified_rank
        if check.get("target_Eh") == target and check.get("fits"):
            uncertainty[j] = check["fits"][j]["estimated_uncertainty_Eh"]*1e6
    levels = {}
    for file in point.glob(sector+"-*.json"):
        if file.name.endswith("-fit.json"):
            continue
        stem = file.stem[len(sector)+1:].removesuffix("-pt2")
        try:
            eps = float(stem)
        except ValueError:
            continue
        if eps not in levels or file.stem.endswith("-pt2"):
            levels[eps] = file
    tightest = min(levels) if levels else None
    raw = read_json(levels[tightest], {}) if levels else {}
    active = [j for j in jobs if j["space"]=="full" and j["r"]==distance and j["sector"]==sector]
    failure = (distance,sector) in failures and not active and not all(passed)
    status = ("done" if all(passed) else "failed" if failure else "running" if active
              else "rank check" if verified_rank is False else "partial" if any(passed)
              else "refining" if current else "sampled" if levels else "queued")
    return dict(status=status,passed=passed,u=uncertainty,levels=len(levels),eps=tightest,
                ndet=raw.get("ndet"),fit_eps=fit.get("epsilon_var") if current else None,
                active=active,rank_check=verified_rank,buffer_eps=check.get("epsilon_var") if check else None)


def collect(root, host, target, live):
    jobs = jobs_on_host(root) if live else []
    failures = set()
    if live:
        name = "hbci-strict-original.log" if host==HOSTS[0] else "hbci-strict-extension.log"
        for line in tail(Path.home()/"MBRefData-runs"/name, 16000):
            match = re.match(r"Failed \(([\d.]+), '([^']+)'\)",line)
            if match:
                failures.add((float(match[1]),match[2]))
    points = []
    for r in host_grid(root,host):
        checks = rank_checks(root,r)
        blocks = [sector_status(root,r,s,target,jobs,checks.get(s),failures) for s in SECTORS]
        points.append(dict(r=r,host=host,sectors=blocks))
    available = None
    if live:
        mem = Path("/proc/meminfo").read_text()
        available = round(int(re.search(r"MemAvailable:\s+(\d+)",mem)[1])/1024**2)
    controls = read_json(root/"calibration/results.json",{})
    validation = read_json(root/"validation/results.json",{})
    return dict(host=host,points=points,jobs=jobs,when=time.time(),live=live,
                load=round(os.getloadavg()[0],1) if live else None,free=available,
                calibration=sum(len(b["fits"]) for b in controls.get("controls",[]) if b.get("qualified"))
                if controls.get("target_Eh")==target else None,
                validation=sum(bool(s["within_target"]) for s in validation.get("states",[]))
                if validation.get("target_Eh")==target else None)


def remote(host,target):
    cmd=["ssh","-o","BatchMode=yes","-o","ConnectTimeout=7","-o","ServerAliveInterval=5",
         "-o","ServerAliveCountMax=1",host,"python3","-","--collect",host,"--target",str(target)]
    result=subprocess.run(cmd,input=SCRIPT.read_bytes(),capture_output=True,timeout=25)
    if result.returncode:
        lines=result.stderr.decode(errors="replace").strip().splitlines()
        raise RuntimeError(lines[-1] if lines else f"SSH exited {result.returncode}")
    data=json.loads(result.stdout)
    if data.get("host")!=host:
        raise ValueError("Unexpected host response")
    return data


def number(x):
    if x is None: return "—"
    if x < .001: return "<.001"
    if x < 1: return f"{x:.3f}".lstrip("0")
    if x < 10: return f"{x:.2f}"
    return f"{x:.0f}"


def epsilon(x):
    return "—" if x is None else f"{x:.1e}".replace("e-0","e-")


def elapsed(sec):
    return f"{int(sec)//3600}h{int(sec)%3600//60:02d}m" if sec>=3600 else f"{int(sec)//60}m"


def table_row(p,wide):
    pieces=[]
    for s in p["sectors"]:
        status=s["status"]
        eps=min((j["eps"] for j in s["active"]),default=s["eps"])
        u="/".join(number(x)+("?" if i==1 and s["rank_check"] is False else "") for i,x in enumerate(s["u"]))
        pieces.append(f"{status:<10} {epsilon(eps):>7} {u:>13}" if wide else f"{status:<10} {u:>13}")
    return f"{p['r']:>4.2f}  {p['host'][-4:]}   "+"  │  ".join(pieces)


def header_row(wide):
    return ("R (Å) host   ¹Ag status    εvar          u₁/u₂  │  ³B₁ᵤ status  εvar          u₁/u₂" if wide else
            "R (Å) host   ¹Ag status            u₁/u₂  │  ³B₁ᵤ status          u₁/u₂")


def host_line(d):
    if d.get("error"):
        kind="last live read" if d["live"] else "local checkout"
        return f"{d['host']}  OFFLINE · {kind} · {d['error']}"
    if not d["live"]:
        return f"{d['host']}  local checkout"
    return (f"{d['host']}  {len(d['jobs'])} active solvers · load {d['load']:.1f} · "
            f"{d['free']} GiB available · read {time.strftime('%H:%M:%S',time.localtime(d['when']))}")


def details(p):
    rows=[]
    for label,b in zip(LABELS,p["sectors"]):
        count=f"{b['ndet']:,}" if b["ndet"] is not None else "—"
        verdict="  rank 2 awaiting independent check" if b["rank_check"] is False else ""
        rows.append(f"{label}  {sum(b['passed'])}/2 qualified · {b['levels']} levels · {count} determinants · fit ε {epsilon(b['fit_eps'])}{verdict}")
        for j in b["active"]:
            tag=j["kind"]
            rows.append(f"  {tag} · ε {epsilon(j['eps'])} · {elapsed(j['elapsed'])} · {j['stage']}")
    return rows


class Monitor:
    def __init__(self,root,published):
        self.root=root
        summary=read_json(root/"results.json")
        self.target=summary["target_per_root_Eh"]
        self.total=summary["requested_roots"]
        self.published=published
        self.pool=concurrent.futures.ThreadPoolExecutor(max_workers=2)
        self.data={h:collect(root,h,self.target,False) for h in HOSTS}
        self.futures={}

    def start(self):
        if self.futures: return
        for h in HOSTS:
            self.futures[h]=(self.pool.submit(collect,self.root,h,self.target,False) if self.published
                             else self.pool.submit(remote,h,self.target))

    def poll(self):
        for h,f in list(self.futures.items()):
            if f.done():
                try: self.data[h]=f.result()
                except (OSError,ValueError,RuntimeError,subprocess.SubprocessError) as exc:
                    self.data[h]["error"]=str(exc)
                del self.futures[h]

    def points(self):
        return sorted((p for d in self.data.values() for p in d["points"]),key=lambda p:p["r"])

    def headline(self):
        n=sum(sum(s["passed"]) for p in self.points() for s in p["sectors"])
        return f"N₂  ·  6-31G HBCI     {n}/{self.total} qualified     target {self.target*1e6:.2f} µEh"

    def controls(self):
        a=self.data[HOSTS[0]]["calibration"]
        b=self.data[HOSTS[1]]["validation"]
        return f"Exact controls   calibration {a if a is not None else '—'}/8    4.00 Å validation {b if b is not None else '—'}/4"

    def close(self):
        self.pool.shutdown(wait=False)


def snapshot(m):
    m.start()
    concurrent.futures.wait(m.futures.values())
    m.poll()
    print(m.headline())
    for h in HOSTS: print(host_line(m.data[h]))
    print(m.controls(),"\n")
    print(header_row(True),"   uncertainties in µEh")
    for p in m.points(): print(table_row(p,True))
    print("\nActive calculations")
    for d in m.data.values():
        for j in d["jobs"]:
            tag=j["kind"]
            print(f"{d['host']}  {j['r']:.2f} Å  {j['sector']}  {tag}  ε {epsilon(j['eps'])}  {elapsed(j['elapsed'])}  {j['stage']}")


def interactive(stdscr,m,interval):
    try: curses.curs_set(0)
    except curses.error: pass
    stdscr.timeout(350)
    stdscr.keypad(True)
    colors={}
    if curses.has_colors():
        curses.start_color()
        try: curses.use_default_colors()
        except curses.error: pass
        for n,(key,col) in enumerate([("title",curses.COLOR_CYAN),("ok",curses.COLOR_GREEN),
                                      ("warn",curses.COLOR_YELLOW),("bad",curses.COLOR_RED)],1):
            try:
                curses.init_pair(n,col,-1)
                colors[key]=curses.color_pair(n)
            except curses.error: pass
    pos=0; top=0; active_only=False; expanded=False; autofocus=True; deadline=0
    while True:
        now=time.monotonic()
        if now>=deadline and not m.futures:
            m.start(); deadline=now+interval
        m.poll()
        pts=m.points()
        if active_only: pts=[p for p in pts if any(b["active"] for b in p["sectors"])]
        if autofocus and not m.futures and not m.published:
            pos=next((i for i,p in enumerate(pts) if any(b["active"] for b in p["sectors"])),pos)
            autofocus=False
        pos=max(0,min(pos,len(pts)-1))
        h,w=stdscr.getmaxyx()
        stdscr.erase()
        def put(y,msg,attr=0):
            if 0<=y<h and w>2:
                try: stdscr.addnstr(y,0,msg,w-1,attr)
                except curses.error: pass
        put(0,m.headline(),curses.A_BOLD|colors.get("title",0))
        for i,host in enumerate(HOSTS):
            d=m.data[host]
            put(i+1,host_line(d),colors.get("warn",0) if d.get("error") else 0)
        put(3,m.controls())
        wide=w>=86
        put(5,header_row(wide),curses.A_BOLD)
        height=1 if expanded else max(1,h-14)
        top=max(0,min(top,max(0,len(pts)-height)))
        if pos<top: top=pos
        if pos>=top+height: top=pos-height+1
        for j,p in enumerate(pts[top:top+height]):
            attr=curses.A_REVERSE if j+top==pos else 0
            if all(all(b["passed"]) for b in p["sectors"]): attr|=colors.get("ok",0)
            elif any(b["status"]=="failed" for b in p["sectors"]): attr|=colors.get("bad",0)
            put(6+j,table_row(p,wide),attr)
        if pts:
            p=pts[pos]
            y=6+min(len(pts),height)+1
            put(y,f"{p['r']:.2f} Å  ·  {p['host']}  ·  uncertainties in µEh",curses.A_BOLD|colors.get("title",0))
            lines=details(p)
            if expanded:
                for j in m.data[p["host"]]["jobs"]:
                    if j["r"]==p["r"]:
                        lines.extend([f"PID {j['pid']}  {j['kind']}  {j['sector']}  {elapsed(j['elapsed'])}",
                                      "  "+j["stage"],"  "+j["log"].split("/MBRefData-runs/")[-1]])
            for i,line in enumerate(lines,1):
                if y+i<h-1: put(y+i,line)
        else:
            put(6,"No active solvers visible in this view.")
        mode="show all" if active_only else "filter active"
        refresh="reading…" if m.futures else f"refresh {max(0,int(deadline-now))}s"
        put(h-1,f"↑/↓ select  Enter details  a {mode}  r refresh  q quit   {refresh}",curses.A_DIM)
        stdscr.refresh()
        key=stdscr.getch()
        if key!=-1: autofocus=False
        if key in (ord("q"),ord("Q"),27): break
        if key in (10,13,curses.KEY_ENTER): expanded=not expanded
        if key in (curses.KEY_DOWN,ord("j")): pos+=1
        elif key in (curses.KEY_UP,ord("k")): pos-=1
        elif key==curses.KEY_NPAGE: pos+=height
        elif key==curses.KEY_PPAGE: pos-=height
        elif key==curses.KEY_HOME: pos=0
        elif key==curses.KEY_END: pos=len(pts)-1
        elif key in (ord("a"),ord("A")): active_only=not active_only; pos=top=0
        elif key in (ord("r"),ord("R")): deadline=0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--once",action="store_true",help="Print one status report and exit")
    parser.add_argument("--published",action="store_true",help="Show the local data checkout without SSH")
    parser.add_argument("--interval",type=int,default=30,help="Seconds between refreshes (default 30)")
    parser.add_argument("--collect",choices=HOSTS,help=argparse.SUPPRESS)
    parser.add_argument("--target",type=float,help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.collect:
        print(json.dumps(collect(Path.home()/"MBRefData-runs"/DATASET,args.collect,args.target,True)))
        return
    if args.interval<5: parser.error("Choose a refresh interval of at least 5 seconds")
    m=Monitor(SCRIPT.parent,args.published)
    try:
        if args.once or not sys.stdout.isatty(): snapshot(m)
        else: curses.wrapper(interactive,m,args.interval)
    except KeyboardInterrupt: pass
    finally: m.close()


if __name__=="__main__":
    main()
