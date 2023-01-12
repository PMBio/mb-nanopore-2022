import tqdm
import gzip
from pathlib import Path
from mb_analysis.config import module_config

from multiprocessing import Process, Queue

def work_readdb(queue, reads, infiles):
    for fq_in in infiles:
        for line in open(str(fq_in) + ".index.readdb", "r"):
            read = line.split("\t")[0]
            if read in reads:
                queue.put(line)
                reads.remove(read)
                if len(reads) == 0:
                    break
        if len(reads) == 0:
            break
    queue.put(None)


def work_fq(queue, reads, infiles):
    for fq_in in infiles:
        with open(fq_in, "r") as f:
            while True:
                fq_lines = [f.readline() for _ in range(4)]
                if fq_lines[0] =="":
                    break
                
                read = fq_lines[0].strip().split(" ")[0][1:]
                if read in reads:
                    queue.put(fq_lines)
                    reads.remove(read)
    queue.put(None)

def main():
    readnames_file = f"{module_config.mbdir}double_minute_2parts/double_minute_readnames.txt"
    fq_in_dir = f"{module_config.mbdir}fastq/Primary/"
    fq_out_file = f"{module_config.mbdir}double_minute_2parts/fastq/dm/0.fq"
    fq_out_readdb_file = f"{module_config.mbdir}double_minute_2parts/fastq/dm/0.fq.index.readdb"
    
    reads = set(r.strip() for r in open(readnames_file, "r"))

    infiles = [fq_in for fq_in in Path(fq_in_dir).iterdir() if fq_in.name.endswith(".fq")]
    
    outqueue = Queue()
    np = 10
    ps = [Process(target=work_readdb, args=(outqueue, reads, infiles[i::np])) for i in range(np)]
    for p in ps:
        p.start()
    
    with open(fq_out_readdb_file, "w", buffering=100) as fq_out_readdb, tqdm.tqdm(total=len(reads)) as pbar:
        while np > 0:
            out_line = outqueue.get()
            
            if out_line is None:
                np -= 1
                continue

            fq_out_readdb.write(out_line)
            pbar.update(1)

    for p in ps:
        p.join()

    np = 10
    ps = [Process(target=work_fq, args=(outqueue, reads, infiles[i::np])) for i in range(np)]
    for p in ps:
        p.start()

    with open(fq_out_file, "w", buffering=100) as fq_out, tqdm.tqdm(total=len(reads)) as pbar:
        while np > 0:
            fqlines = outqueue.get()
            if fqlines is None:
                np -= 1
                print("One process finished. Left: ", np)
                continue
            for outline in fqlines:
                fq_out.write(outline)
            pbar.update(1)
       
    for p in ps:
        p.join()

    with open(fq_out_file, "w") as fq_out, tqdm.tqdm(total=len(reads)) as pbar:
        for fq_in in infiles:
            with open(fq_in, "r") as f:
                while True:
                    fq_lines = [f.readline() for _ in range(4)]
                    if fq_lines[0] =="":
                        break
                    print(fq_lines[0])
                    read = fq_lines[0].strip().split(" ")[0][1:]
                    if read in reads:
                        for outline in fq_lines:
                            fq_out.write(outline)
                        reads.remove(read)
                        pbar.update(1)
                        if len(reads) == 0:
                            break
                if len(reads) == 0:
                    break

if __name__ == "__main__":
    main()