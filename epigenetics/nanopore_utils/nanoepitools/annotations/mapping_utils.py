from nanoepitools.annotations.annotations import GFFFeature, GFFAnnotationsReader


class MapToPromoter:
    """
    Maps a GFF feature to a promoter it overlaps with
    """
    def __init__(self, gff: GFFAnnotationsReader, promoter_before_tss: int, promoter_after_tss: int, input_will_be_sorted=False):
        self.gff = gff
        self.promoter_before_tss = promoter_before_tss
        self.promoter_after_tss = promoter_after_tss
        self.input_will_be_sorted = input_will_be_sorted
        self.last_chrom = None
        self.seeker = None
    
    def range_transform(self, gff_feature: GFFFeature):
        if gff_feature.direction == "+":
            promoter_range = [gff_feature.start - self.promoter_before_tss, gff_feature.start + self.promoter_after_tss]
        elif gff_feature.direction == "-":
            promoter_range = [gff_feature.end - self.promoter_after_tss, gff_feature.end + self.promoter_before_tss]
        
        return promoter_range
    
    def dist_function(self, gff_feature: GFFFeature, start: int, end: int):
        """Computes the distance of the region to the promoter of the given gff feature.
        Returns 0 as a distance if there is overlap (i.e. the promoter is fully or partially contained within the
        region or vice versa)
        """
        
        """Promoter location depends on the direction of the transcript"""
        promoter_range = self.range_transform(gff_feature)
        
        if promoter_range[0] > end or promoter_range[1] < start:
            """Our region is close to the promoter but not overlapping"""
            return 1
        else:
            """Our region is overlapping the promoter, therefore the distance is 0"""
            return 0
    
    def find_sorted(self, chrom, start, end):
        if chrom != self.last_chrom and self.input_will_be_sorted:
            self.last_chrom = chrom
            self.seeker = self.gff.chromosomes[chrom].get_sorted_in_range_finder()
        return self.seeker.find(start, end, range_transform=lambda x: self.range_transform(x), max_recursion=1, min_recursion=1)
    
    def find_unsorted(self, chrom, start, end):
        return self.gff.chromosomes[chrom].get_in_range(start, end, range_transform=lambda x: self.range_transform(x),
            max_recursion=1, min_recursion=1)
    
    def __call__(self, chrom, start, end):
        if self.input_will_be_sorted:
            return self.find_sorted(chrom, start, end)
        else:
            return self.find_unsorted(chrom, start, end)
