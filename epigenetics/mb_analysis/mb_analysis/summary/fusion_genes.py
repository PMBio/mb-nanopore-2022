from typing import List
import re

import pandas as pd
import mygene
import itertools

from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature


def has_gene_id(geneidlist: List[str]):
    return any([geneid is not None for geneid in geneidlist])


class FusionPartner:
    def __init__(self, gene, confidence):
        self.gene = gene
        self.confidence = confidence
    
    def __repr__(self):
        return f"{self.gene} ({self.confidence})"
    
    def __str__(self):
        return self.__repr__()
    
    def __hash__(self):
        return self.gene.__hash__()
    
    def __eq__(self, other):
        return self.gene == other.gene


class FusionGenes:
    """
    This class loads the Arriba output and has some helper functions for mapping gene names
    using the mygene python package
    """
    def __init__(self, gff: GFFAnnotationsReader, verbose=False):
        self.fusion_df_unfiltered = None
        self.gene_name_to_id_lookup = {
            gene.name: gene.sanitized_id() for chrom in gff.chromosomes.values() for gene in chrom.children.values()
        }
        self.fusion_pairs = {}
        self.genename_regex = re.compile("^([A-Za-z0-9-\.]*)[/(]*")
        self.mg = mygene.MyGeneInfo()
        self.verbose = verbose
    
    def __lookup_gene_id_from_gff(self, genename: str) -> str:
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
        
        id = self.__lookup_gene_id_from_gff(r["query"])
        if id is not None:
            return id
        # It wasn't in the gff, check the mg results
        
        if "ensembl" in r:
            # we have results from the original query
            return self.__resolve_query_result(r["ensembl"])
        else:
            # It wasn't in the GFF, try an individual, broader search
            new_r = self.mg.query(r["query"], fields="ensembl.gene", species="human")["hits"]
            if len(new_r) == 0:
                if self.verbose:
                    print("No result", r["query"])
                return None
            for new_r_hit in new_r:
                if "ensembl" not in new_r_hit:
                    continue
                id = self.__resolve_query_result(new_r_hit["ensembl"])
                if id is not None:
                    return id
            if self.verbose:
                print("No result", r["query"])
            return None
    
    def __lookup_gene_ids(self, genesymbols):
        symbol_cleaned_dict = {
            symbol: self.genename_regex.match(symbol).group(1)
            for row_symbols in genesymbols
            for symbol in row_symbols.split(",")
        }
        cleaned_symbol_dict = {v: k for k, v in symbol_cleaned_dict.items()}
        query_result = self.mg.querymany(
            list(symbol_cleaned_dict.values()), fields="ensembl.gene", species="human", scopes="symbol"
        )
        mapping_dict = {}
        for r in query_result:
            table_symbol = cleaned_symbol_dict[r["query"]]
            mapping_dict[table_symbol] = self.__research_if_not_found(r)
        return genesymbols.map(
            lambda row_symbols: [mapping_dict.get(symbol, None) for symbol in row_symbols.split(",")]
        )
    
    def load(self, filepath: str):
        self.fusion_df = pd.read_csv(filepath, sep="\t")
        for column in self.fusion_df.columns:
            if set(self.fusion_df[column]) == {0, 1}:
                self.fusion_df[column] = self.fusion_df[column].astype(bool)
        
        self.fusion_df["geneid1"] = self.__lookup_gene_ids(self.fusion_df["X.gene1"])
        self.fusion_df["geneid2"] = self.__lookup_gene_ids(self.fusion_df["gene2"])
        self.fusion_df_unfiltered = self.fusion_df
        return self
    
    def filter(self, column=None, mode=None):
        if column is not None:
            self.fusion_df = self.fusion_df.loc[self.fusion_df[column] == 1]
        if mode == "has_gene_id":
            self.fusion_df = self.fusion_df.loc[self.fusion_df["geneid1"].map(has_gene_id)]
            self.fusion_df = self.fusion_df.loc[self.fusion_df["geneid2"].map(has_gene_id)]
        if mode == "high_confidence":
            self.fusion_df = self.fusion_df.loc[self.fusion_df["confidence"] == "high"]
        return self
    
    def reset_filter(self):
        self.fusion_df = self.fusion_df_unfiltered
    
    def __iter__(self):
        return (
            (perm_a, FusionPartner(perm_b, row["confidence"]))
            for _, row in self.fusion_df.iterrows()
            for gene_a in row["geneid1"]
            for gene_b in row["geneid2"]
            for perm_a, perm_b in itertools.permutations([gene_a, gene_b])
            if perm_a is not None and perm_b is not None
        )
