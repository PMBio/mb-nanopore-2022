from pathlib import Path

CHAINFILE_COL_QUERY_CHROM = 2
CHAINFILE_COL_QUERY_DIRECTION = 4
CHAINFILE_COL_QUERY_START = 5
CHAINFILE_COL_QUERY_END = 6

CHAINFILE_COL_TARGET_CHROM = 7
CHAINFILE_COL_TARGET_DIRECTION = 9
CHAINFILE_COL_TARGET_START = 10
CHAINFILE_COL_TARGET_END = 11


class Chain:
    def __init__(
        self, ref_contig, ref_start, ref_end, ref_direction, query_contig, query_start, query_end, query_direction
    ):
        self.ref_contig = ref_contig
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.ref_direction = ref_direction
        self.query_contig = query_contig
        self.query_start = query_start
        self.query_end = query_end
        self.query_direction = query_direction
        self.lines = []

    def __iter__(self):
        return iter(self.lines)

    def __str__(self):
        return f"Ref: {self.ref_contig}:{self.ref_start}-{self.ref_end} ({self.ref_direction}) Query: {self.query_contig}:{self.query_start}-{self.query_end} ({self.query_direction})"

    def translate_sorted_target_to_query(self, pos_list):
        query_offset = 0
        ref_offset = 0
        result = []
        if self.ref_direction == "-":
            print("rev")
            # Reading in reverse
            pos_list = [self.ref_end + self.ref_start - p for p in pos_list][::-1]
        pos_list = iter(pos_list)
        pos = next(pos_list)
    
        try:
            for link in self:
                while True:  # Until no position matched here
                    if pos < self.ref_start + ref_offset:
                        # Next match-group starts after requested position
                        result.append(None)
                        pos = next(pos_list)
                        continue
                    ref_end_match = self.ref_start + ref_offset + link[0]
                    diff_start = ref_end_match - pos
                    # test if pos < ref_end_match
                    if diff_start > 0:
                        # Requested position is in this match-group
                        result.append((self.query_contig, self.query_start + query_offset + link[0] - diff_start))
                        pos = next(pos_list)
                        continue
                    if link[1] is not None:
                        # If this isnt the last link, prepare to look at next link
                        ref_offset += link[2] + link[0]
                        query_offset += link[1] + link[0]
                    break
        except StopIteration:
            pass
        if self.ref_direction == "-":
            return result[::-1]
        else:
            return result

class ChainSet:
    def __init__(self):
        self.chains = []

    def __iter__(self):
        return iter(self.chains)

    def translate_target_to_query(self, pos):
        for chain in self:
            if chain.ref_start < pos < chain.ref_end:
                query_offset = 0
                ref_offset = 0
                for link in chain:
                    if chain.ref_start + ref_offset > pos:
                        return None
                    ref_end_match = chain.ref_start + ref_offset + link[0]
                    diff_start = ref_end_match - pos
                    if diff_start >= 0:
                        return (chain.query_contig, chain.query_start + query_offset + link[0] - diff_start)
                    else:
                        if link[1] is None:
                            return None
                        ref_offset += link[1] + link[0]
                        query_offset += link[2] + link[0]

    def translate_sorted_target_to_query(self, pos_list):
        pos_enum = enumerate(pos_list)
        i, pos = next(pos_enum)
    
        chain_iter = iter(self.chains)
        chain = next(chain_iter)
    
        combined_mapping = []
    
        while True:
            if chain.ref_start <= pos < chain.ref_end:
                # Find remaining positions:
                i_start = i
                while True:
                    if i == len(pos_list) - 1:
                        # Last position
                        i = i + 1
                        break
                    i, pos = next(pos_enum)
                    if pos > chain.ref_end:
                        break
                mapping = chain.translate_sorted_target_to_query(pos_list[i_start: i])
                combined_mapping = combined_mapping + mapping
        
            if i == len(pos_list):
                # Last position
                break
        
            if pos < chain.ref_start:
                combined_mapping.append(None)
                i, pos = next(pos_enum)
                continue
            if chain.ref_end <= pos:
                chain = next(chain_iter)
                continue
        return combined_mapping



class ChainReader:
    def __init__(
        self,
        parent,
        ref_contig,
        ref_start,
        ref_end,
        ref_direction,
        query_contig,
        query_start,
        query_end,
        query_direction,
    ):
        self.parent = parent
        self.ref_contig = ref_contig
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.ref_direction = ref_direction
        self.query_contig = query_contig
        self.query_start = query_start
        self.query_end = query_end
        self.query_direction = query_direction
        self.done = False
        self.lines = None

    def __iter__(self):
        return self

    def read_next(self):
        if self.done:
            # Makes sure we don't read over the chain
            raise StopIteration

        line = next(self.parent.f).strip()

        if len(line) == 0:
            self.done = True
            raise StopIteration

        ret = [int(r) for r in line.split("\t")]

        if len(ret) == 3:
            return ret

        elif len(ret) == 1:
            return [ret[0], None, None]

    def __next__(self):
        return self.read_next()


class ChainFileReader:
    def __init__(self, filepath: Path):
        self.f = open(filepath, "r")
        self.current_chain_reader = None

    def __enter__(self, *argc):
        return self

    def __iter__(self):
        return self

    def __exit__(self, *argc):
        self.f.__exit__(*argc)

    def __next__(self):
        if self.current_chain_reader is not None:
            # Make sure we are done reading the last chain
            for _ in self.current_chain_reader:
                pass

        while True:
            line = next(self.f).strip()
            if not line.startswith("#"):
                break

        assert line.startswith("chain")
        line = line.split(" ")
        chainreader = ChainReader(
            self,
            line[CHAINFILE_COL_TARGET_CHROM],
            int(line[CHAINFILE_COL_TARGET_START]),
            int(line[CHAINFILE_COL_TARGET_END]),
            line[CHAINFILE_COL_TARGET_DIRECTION],
            line[CHAINFILE_COL_QUERY_CHROM],
            int(line[CHAINFILE_COL_QUERY_START]),
            int(line[CHAINFILE_COL_QUERY_END]),
            line[CHAINFILE_COL_QUERY_DIRECTION],
        )
        self.current_chain_reader = chainreader
        return chainreader

    def read_to_memory(self):
        chain_set = ChainSet()
        for chain_reader in self:
            chain = Chain(
                chain_reader.ref_contig,
                chain_reader.ref_start,
                chain_reader.ref_end,
                chain_reader.ref_direction,
                chain_reader.query_contig,
                chain_reader.query_start,
                chain_reader.query_end,
                chain_reader.query_direction,
            )
            chain.lines = [line for line in chain_reader]
            chain_set.chains.append(chain)
        return chain_set
