from __future__ import annotations
import re
from typing import Dict, List
import itertools

import mygene
import numpy as np


def range_comparison_dist_fn(feature: GFFFeature, start, end):
    """Distance function between """
    if feature.start > end or feature.end < start:
        # Not overlapping:
        return min(abs(feature.end - start), abs(feature.start - end))
    else:
        # Overlapping:
        return 0


class FeatureInRangeFinder:
    """This class shouldn't be instantiated outside of the GFFFeatures
    It is used to avoid copy&pasting the code for finding features that overlap
    a certain genomic range in (a) sorted and (b) unsorted lists"""
    
    def __init__(self, parent_feature: GFFFeature, input_is_sorted: bool):
        """
        :param parent_feature: Typically a chromosome feature
        :param input_is_sorted: Whether the ranges that will be queried later will be sorted.
                                If this is true, the search will be O(N+M) time, whereas
                                if this is false, the search will be N*M time. But beware,
                                there are no sanity checks! If you say they are sorted, and
                                then give me unsorted ranges, you will miss hits!
        """
        self.parent_feature = parent_feature
        self.earliest_index = 0
        self.input_is_sorted = input_is_sorted
    
    def find(self, start: int, end: int, range_transform=None, max_recursion=0, min_recursion=0):
        """
        Returns a generator! Finds all child features which overlap [start,end]
        :param start: the start location
        :param end:
        :param range_transform:
        :param max_recursion:
        :param min_recursion:
        :return:
        """
        children = self.parent_feature.sorted_children
        if range_transform is None:
            # identity
            def range_transform(x):
                return x.start, x.end
        
        else:
            # If range transform is not identity, we need to sort these accordingly,
            # or else we can no longer assume they are ordered
            children = sorted(children, key=lambda x: range_transform(x)[0])
        if self.input_is_sorted:
            children = children[self.earliest_index :]
        for feature in children:
            if range_transform(feature)[0] > end and min_recursion <= 0:
                break
            elif range_transform(feature)[1] < start and min_recursion <= 0:
                if self.input_is_sorted:
                    self.earliest_index += 1
                continue
            else:
                if max_recursion > 0:
                    for x in feature.in_range_finder.find(
                        start,
                        end,
                        range_transform=range_transform,
                        max_recursion=max_recursion - 1,
                        min_recursion=min_recursion - 1,
                    ):
                        yield x
                else:
                    yield feature


class GFFFeature:
    def __init__(self, start, end, type, id, direction, name=None, parent=None):
        self.start = start
        self.end = end
        self.type: str = type
        self.id: str = id
        self.direction = direction
        if name is None and id is not None:
            if ":" in id:
                name = id.split(":")[1]
            else:
                name = id
        self.name: str = name
        self.children: Dict[str, GFFFeature] = {}
        self.leafs: List[GFFFeature] = []
        
        self.sorted_children: List[GFFFeature] = []
        self.parent = parent
        # Default finder
        self.in_range_finder = FeatureInRangeFinder(self, input_is_sorted=False)
    
    def sanitized_id(self):
        if ":" in self.id:
            return self.id.split(":")[1]
        else:
            return self.id
    
    def build_sorted(self):
        self.sorted_children = sorted(itertools.chain(self.children.values(), self.leafs), key=lambda x: x.start)
        for child in self.sorted_children:
            child.build_sorted()
    
    def get_in_range(self, start, end, range_transform=None, min_recursion=0, max_recursion=0):
        return self.in_range_finder.find(
            start, end, range_transform=range_transform, min_recursion=min_recursion, max_recursion=max_recursion
        )
    
    def get_sorted_in_range_finder(self):
        """If you are searching for annotations for a SORTED list of ranges,
        then this allows you to do so in a near linear time, as it continues
        to move up the pointer at which searching starts"""
        
        return FeatureInRangeFinder(self, input_is_sorted=True)
    
    def get_nearest_feature(
        self, start=None, end=None, dist_fn=None, nearest=None, nearest_dist=10e15, max_recursion=0, min_recursion=0
    ):
        if dist_fn is None:
            dist_fn = lambda x: range_comparison_dist_fn(x, start, end)
        
        for feature in self.sorted_children:
            if min_recursion > 0 and len(feature.sorted_children) > 0:
                _, dist = feature.get_nearest_feature(
                    dist_fn=dist_fn, max_recursion=max_recursion - 1, min_recursion=min_recursion - 1
                )
            else:
                dist = dist_fn(feature)
            
            if max_recursion > 0 and len(feature.sorted_children) == 0:
                # If we are looking for a deeper feature but this one doesn't contain any, skip it
                continue
            
            if dist < nearest_dist:
                nearest_dist = dist
                nearest = feature
            else:
                continue
        
        if nearest is None:
            return None, np.inf
        
        if max_recursion > 0:
            return nearest.get_nearest_feature(
                dist_fn=dist_fn, max_recursion=max_recursion - 1, min_recursion=min_recursion - 1
            )
        else:
            return nearest, nearest_dist
    
    def __repr__(self):
        parent = self
        while parent.parent is not None:
            parent = parent.parent
        chrom = parent.name
        return f"{self.id} ({chrom}:{self.start}-{self.end})"


class GFFAnnotationsReader:
    def __init__(self):
        self.chromosomes: Dict[str, GFFFeature] = {}
        self.included_features = ["gene", "mRNA", "exon", "ncRNA_gene", "pseudogene", "lnc_RNA"]
    
    def read(self, gff_file, only_protein_coding=True, chroms=None):
        with open(gff_file, "rt") as f:
            
            rowit = (
                {"chr": l[0], "type": l[2], "start": int(l[3]), "end": int(l[4]), "direction": l[6], "info": l[8]}
                for l in (l.split("\t") for l in f if l[0] != "#")
            )
            
            cur_chrom = None
            for row in rowit:
                if chroms is not None and row["chr"] not in chroms:
                    continue
                if cur_chrom is None or cur_chrom.id != row["chr"]:
                    cur_chrom = GFFFeature(0, 0, type="chrom", direction="+", id=row["chr"])
                    self.chromosomes[row["chr"]] = cur_chrom
                    open_features: Dict[str, GFFFeature] = {}
                
                if row["type"] not in self.included_features:
                    continue
                
                info = [kv.split("=") for kv in row["info"].split(";")]
                info = {k: v for k, v in info}
                id = info.get("ID", None)
                name = info.get("Name", None)
                parent_id = info.get("Parent", None)
                
                if "biotype" in info:
                    if only_protein_coding and not info["biotype"] == "protein_coding":
                        continue
                elif only_protein_coding:
                    # "only_protein_coding" requires a biotype annotation
                    continue
                
                cur_feature = GFFFeature(
                    row["start"], row["end"], type=row["type"], id=id, direction=row["direction"], name=name
                )
                
                if parent_id is not None:
                    if parent_id not in open_features:
                        continue
                    cur_parent = open_features[parent_id]
                else:
                    cur_parent = cur_chrom
                
                cur_feature.parent = cur_parent
                
                if id is None:
                    # Only leafs have no id
                    cur_parent.leafs.append(cur_feature)
                else:
                    cur_parent.children[id] = cur_feature
                    open_features[id] = cur_feature
    
    def build_index(self):
        for chrom_container in self.chromosomes.values():
            chrom_container.build_sorted()
    
    def get_gene(self, gene_id: str) -> GFFFeature:
        """
        :param gene_id: Gene ID without the "gene:" part at the beginning
        """
        gene_id = f"gene:{gene_id}"
        genes_found = [chrom.children[gene_id] for chrom in self.chromosomes.values() if gene_id in chrom.children]
        assert len(genes_found) <= 1, "Multiple genes with the same id were found"
        if len(genes_found) == 0:
            return None
        else:
            return genes_found[0]
    
    def get_transcript(self, transcript_id: str) -> GFFFeature:
        """
        :param transcript_id: Transcript ID without the "transcript:" part at the beginning
        """
        transcript_id = f"transcript:{transcript_id}"
        for chrom in self.chromosomes.values():
            for gene in chrom.children.values():
                if transcript_id in gene.children:
                    return gene.children[transcript_id]
        return None
    
    def _print(self, feature, depth=2):
        for fid, f in feature.children.items():
            print("".join(["-"] * depth), fid)
            self._print(f, depth + 2)
        for f in feature.leafs:
            print("".join(["-"] * depth), f.start, f.end)
    
    def print(self):
        for chromname, chrom in self.chromosomes.items():
            print(chromname)
            self._print(chrom)


class GeneNameToEnsemblID:
    def __init__(self, gff: GFFAnnotationsReader, verbose=False):
        self.gene_name_to_id_lookup = {
            gene.name: gene.sanitized_id() for chrom in gff.chromosomes.values() for gene in chrom.children.values()
        }
        self.genename_regex = re.compile("^([A-Za-z0-9-\.]*)[/(]*")
        self.mg = mygene.MyGeneInfo()
        self.verbose = verbose
    
    def lookup_gene_id_from_gff(self, genename: str) -> str:
        genename = self.genename_regex.match(genename).group(1)
        return self.gene_name_to_id_lookup.get(genename, None)
    
    def __resolve_query_result(self, r):
        if isinstance(r, list):
            # There was multiple results
            for r_entry in r:
                if r_entry["gene"] in self.gene_name_to_id_lookup.values():
                    # Prefer result that is also in the gff
                    return r_entry["gene"]
            return r[0]["gene"]  # first result as default
        else:
            # There was one result
            return r["gene"]
    
    def __research_if_not_found(self, r):
        if "ensembl" in r:
            # we have results from the original query
            return self.__resolve_query_result(r["ensembl"])
        else:
            return None
    
    def lookup_gene_ids_online(self, genesymbols):
        symbol_cleaned_dict = {
            symbol: self.genename_regex.match(symbol).group(1)
            for row_symbols in genesymbols
            for symbol in row_symbols.split(",")
        }
        cleaned_symbol_dict = {v: k for k, v in symbol_cleaned_dict.items()}
        query_result = self.mg.querymany(
            list(symbol_cleaned_dict.values()),
            fields="ensembl.gene",
            species="human",
            scopes=["symbol", "alias"],
            verbose=self.verbose,
        )
        mapping_dict = {}
        for r in query_result:
            table_symbol = cleaned_symbol_dict[r["query"]]
            mapping_dict[table_symbol] = self.__research_if_not_found(r)
        return mapping_dict
    
    def translate_gene_names_to_ids(self, gene_symbols):
        translated = {name: self.lookup_gene_id_from_gff(name) for name in gene_symbols}
        untranslated = [name for name, id in translated.items() if id is None]
        translated.update(self.lookup_gene_ids_online(untranslated))
        return [translated.get(name, None) for name in gene_symbols]
