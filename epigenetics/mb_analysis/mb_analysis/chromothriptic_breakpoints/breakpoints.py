from __future__ import annotations
import heapq
from typing import List

import numpy as np
import pandas as pd

from mb_analysis.config import module_config

from nanoepitools.annotations.annotations import GFFFeature, GFFAnnotationsReader


def build_nodes(row):
    return {"chrom": row["query.chr"], "start": row["query.start"], "dir": "-" if row["query.ct"][0] == "3" else "+"}, {
        "chrom": row["query.chr2"],
        "start": row["query.end"],
        "dir": "-" if row["query.ct"][3] == "3" else "+",
    }


class BPNode:
    def __init__(self, dist, node):
        self.dist = dist
        self.node = node
    
    def __lt__(self, other):
        return self.dist < other.dist
    
    def __repr__(self):
        return f"(d={self.dist}: node={self.node})"


class ShortestPathResult:
    def __init__(self, is_connected: bool, uses_breakpoint: bool, distance: float, path: List[BPNode]):
        self.is_connected = is_connected
        self.uses_breakpoint = uses_breakpoint
        self.distance = distance
        self.path = path
    
    def __repr__(self):
        return f"ShortestPathResult(is_connected={self.is_connected}, uses_breakpoint={self.uses_breakpoint}, distance={self.distance})"


class ChromothripticBreakpoints:
    """
    This class loads chromothriptic SVs (inter and intra chromosomal) and implements a shortest path algorithm to
    explain how gene fusions could have been created
    """
    def __init__(self, gff: GFFAnnotationsReader):
        self.breakpoint_df: pd.Dataframe = None
        self.bp_nodes = []
        self.bp_edges = {}
        self.gff = gff
        # Gene ids that are allowed to be crossed when computing shortest path (beginning and end node)
        # This is set by get_shortest_path
        self.__crossing_allowed = {}
    
    def __test_gene_crossing(self, chrom, start, end):
        between_exons = {
            g.parent.parent.sanitized_id() for g in self.gff.chromosomes[chrom].get_in_range(start, end, max_recursion=2)
        }
        between_exons = between_exons.difference(self.__crossing_allowed)
        return len(between_exons) > 0
    
    def __node_dist(self, node_a, node_b, max_dist=1e4):
        edge = (frozenset(node_a.items()), frozenset(node_b.items()))
        if edge in self.bp_edges:
            return 1
        
        if node_a["chrom"] != node_b["chrom"]:
            # not the same chrom
            return np.inf
        chrom = node_a["chrom"]
        
        if "end" in node_a and "end" in node_b:
            # both are ranges
            if node_a["end"] < node_b["start"] or node_b["end"] < node_a["start"]:
                # no overlap
                start = min(node_a["end"], node_b["end"])
                end = max(node_a["start"], node_b["start"])
            else:
                # they are overlapping
                return 0
        elif "end" in node_a or "end" in node_b:
            # one of them is a range
            range_node = node_a if "end" in node_a else node_b
            spot_node = node_b if "end" in node_a else node_a
            if range_node["start"] <= spot_node["start"] <= range_node["end"]:
                return 0
            elif spot_node["start"] > range_node["end"]:
                if spot_node.get("dir", "") == "+":
                    # Spot node connects to higher valued range, but range node is before
                    return np.inf
                start = range_node["end"]
                end = spot_node["start"]
            else:
                if spot_node.get("dir", "") == "-":
                    # Spot node connects to lower valued range, but range node is after
                    return np.inf
                start = spot_node["start"]
                end = range_node["start"]
        
        else:
            # neither of them is a range
            if node_a["start"] < node_b["start"]:
                if node_a.get("dir", "+") == "+" and node_b.get("dir", "-") == "-":
                    start = node_a["start"]
                    end = node_b["start"]
                else:
                    return np.inf
            else:
                if node_a.get("dir", "-") == "-" and node_b.get("dir", "+") == "+":
                    start = node_b["start"]
                    end = node_a["start"]
                else:
                    return np.inf
        
        dist = end - start
        if dist > max_dist:
            return np.inf
        
        if self.__test_gene_crossing(chrom, start, end):
            # We are not allowed to cross other genes!
            return np.inf
        return dist
    
    def node_dist(self, *args, **kwargs):
        return self.__node_dist(*args, **kwargs)
    
    def __build_graph(self, max_dist=1e4):
        node_set = set()
        for _, row in self.breakpoint_df.iterrows():
            node_set.update({frozenset(x.items()) for x in build_nodes(row)})
        
        self.bp_nodes = [dict(node) for node in node_set]
        
        for _, row in self.breakpoint_df.iterrows():
            edge = build_nodes(row)
            self.bp_edges[frozenset(edge[0].items()), frozenset(edge[1].items())] = 1
            self.bp_edges[frozenset(edge[1].items()), frozenset(edge[0].items())] = 1
    
    def load(self) -> ChromothripticBreakpoints:
        self.breakpoint_df = pd.read_csv(module_config.ct_breakpoints_file, sep="\t")
        self.breakpoint_df["query.chr"] = self.breakpoint_df["query.chr"].map(lambda x: x.replace("chr", ""))
        self.breakpoint_df["query.chr2"] = self.breakpoint_df["query.chr2"].map(lambda x: x.replace("chr", ""))
        self.__build_graph()
        return self
    
    def __iter__(self):
        return (row for _, row in self.breakpoint_df.iterrows())
    
    def __get_shortest_path(self, start: BPNode, end: BPNode, max_dist=3e4, genes_crossing_allowed=None):
        """
        This method implements Dijkstra's algorithm for shortest path between two nodes in a graph
        """
        if genes_crossing_allowed is not None:
            self.__crossing_allowed = genes_crossing_allowed
        else:
            self.__crossing_allowed = {}
        
        all_nodes = [start] + [BPNode(np.inf, node) for node in self.bp_nodes] + [end]
        q = [start]
        visited = set()
        dist = {node: node.dist for node in all_nodes}
        reach_via = {}
        direct_distance = self.__node_dist(start.node, end.node, max_dist=max_dist)
        
        while len(q) > 0:
            v = heapq.heappop(q)
            if v in visited:
                continue
            visited.update({v})
            
            for z in all_nodes:
                if v == z:
                    continue
                
                d = self.__node_dist(v.node, z.node, max_dist=max_dist)
                
                d = dist[v] + d
                if d > max_dist:
                    # Threshold by maximum distance
                    d = np.inf
                
                if np.isinf(d):
                    # not a connected node
                    continue
                
                if d < dist[z]:
                    dist[z] = d
                    reach_via[z] = v
                    z.dist = d
                    heapq.heappush(q, z)
        
        final_distance = dist[end]
        is_connected = not np.isinf(final_distance)
        shortest_path = []
        if is_connected:
            prev = end
            shortest_path.append(end)
            while prev != start and prev is not None:
                prev = reach_via.get(prev, None)
                shortest_path.append(prev)
        return ShortestPathResult(
            is_connected=is_connected,
            uses_breakpoint=final_distance < direct_distance,
            distance=final_distance,
            path=shortest_path,
        )
    
    def get_shortest_path(self, chr_a: str, pos_a: int, chr_b: str, pos_b: int, genes_crossing_allowed=None, **kwargs):
        start = BPNode(0, {"chrom": chr_a, "start": pos_a})
        end = BPNode(np.inf, {"chrom": chr_b, "start": pos_b})
        
        start_genes = {g.sanitized_id() for g in self.gff.chromosomes[chr_a].get_in_range(pos_a, pos_a + 1)}
        end_genes = {g.sanitized_id() for g in self.gff.chromosomes[chr_b].get_in_range(pos_b, pos_b + 1)}
        
        included_genes = start_genes.union(end_genes)
        if genes_crossing_allowed is None:
            genes_crossing_allowed = included_genes
        else:
            genes_crossing_allowed = genes_crossing_allowed.union(included_genes)
        
        return self.__get_shortest_path(start, end, genes_crossing_allowed=genes_crossing_allowed, **kwargs)
    
    def get_shortest_path_between_two_ranges(
        self,
        chr_a: str,
        start_a: int,
        end_a: int,
        chr_b: str,
        start_b: int,
        end_b: int,
        genes_crossing_allowed=None,
        **kwargs,
    ):
        start_genes = {g.sanitized_id() for g in self.gff.chromosomes[chr_a].get_in_range(start_a, end_a)}
        end_genes = {g.sanitized_id() for g in self.gff.chromosomes[chr_b].get_in_range(start_b, end_b)}
        
        included_genes = start_genes.union(end_genes)
        if genes_crossing_allowed is None:
            genes_crossing_allowed = included_genes
        else:
            genes_crossing_allowed = genes_crossing_allowed.union(included_genes)
        
        start = BPNode(0, {"chrom": chr_a, "start": start_a, "end": end_a})
        end = BPNode(np.inf, {"chrom": chr_b, "start": start_b, "end": end_b})
        return self.__get_shortest_path(start, end, genes_crossing_allowed=genes_crossing_allowed, **kwargs)
    
    def get_shortest_path_between_two_genes(self, gene_a: GFFFeature, gene_b: GFFFeature, **kwargs):
        return self.get_shortest_path_between_two_ranges(
            gene_a.parent.sanitized_id(), gene_a.start, gene_a.end, gene_b.parent.sanitized_id(), gene_b.start, gene_b.end, **kwargs
        )
