import pandas as pd


def intersect_locations_with_ranges(locations, ranges, nonevalue=-1):
    """
    Takes a list of sites, and a list of genomic regions, and maps
    the sites to the regions in linear time (O(n+m), where n is the number
    of locations and m the number of regions). Requires both locations and
    ranges to be sorted.

    A few important prerequisites on the input:
      * Both dataframes need to have the columns in the order
        specified below, since the function uses itertuples to iterate
        quickly. Column names are irrelevant.
      * Both dataframes need to be sorted by the chromosome column, such
        that we can compare the two. So, if they are both using strings
        for the chromosome, they should both be sorted lexographically.
        If you are unsure about this, set is_sorted to false, so this
        function will perform the correct sorting for you.
      * Both dataframes are allowed to have further columns, as long
        as the first few columns are as specified.

    :param locations: pandas dataframe of the format: (location, ...)
    :param ranges: pandas dataframe of the format:
                   (start_location, end_location,...)
    :param nonevalue: what value to use for a "None" index. Since not all
        datatypes are nullable, you should choose whatever works best
        for your datatype. Default is -1 since this works with int64.
        CAVEAT: If you provide "None" but the index datatype is not
        nullable, pandas will just convert it (to float, for example)
        and you will end up with indices that are not comparable.
    :return: A panda series where the index is the same index as
             in locations, and the value is either the index of the matching
             region, or None if the location is in no region
    """

    region_membership = pd.Series(
        data=nonevalue, index=locations.index, dtype=ranges.index.dtype
    )

    en_ranges = ranges.itertuples()
    en_loc = locations.itertuples()

    region = next(en_ranges)
    loc = next(en_loc)

    try:
        """
        When accessing the tuples, remember:
            loc[0]: cpg index
            loc[1]: location
            region[0]: region index
            region[1]: start site
            region[2]: end site
        """
        while True:
            """ If location is behind region, spool location """
            while loc[1] < region[1]:
                loc = next(en_loc)

            """ If location is past region, spool region """
            while loc[1] > region[2]:
                region = next(en_ranges)

            """ Check all the constraints """
            if region[1] <= loc[1] <= region[2]:
                region_membership[loc[0]] = region[0]
            else:
                """ This happens if we spooled the region past the location """
                pass

            """ Get next cpg location """
            loc = next(en_loc)

    except StopIteration:
        pass

    return region_membership


def intersect_locations_with_ranges_by_chromosome(locations, ranges, locations_chrom_key="chrom",
                                                  ranges_chrom_key="chrom", nonevalue=-1):
    chromosome_intersection = set(locations[locations_chrom_key]).intersection(set(ranges[ranges_chrom_key]))
    locations = locations.groupby(locations_chrom_key)
    ranges = ranges.groupby(ranges_chrom_key)
    
    region_membership = []
    for chrom in chromosome_intersection:
        chr_locs = locations.get_group(chrom).drop(locations_chrom_key, axis=1)
        chr_rngs = ranges.get_group(chrom).drop(ranges_chrom_key, axis=1)
        region_membership.append(intersect_locations_with_ranges(chr_locs, chr_rngs, nonevalue=nonevalue))
    
    return pd.concat(region_membership)


def intersect_ranges(ranges_a, ranges_b, nonevalue=-1):
    region_membership = pd.Series(data=nonevalue, index=ranges_a.index, dtype=ranges_b.index.dtype)
    
    en_a = ranges_a.itertuples()
    en_b = ranges_b.itertuples()
    
    a = next(en_a)
    b = next(en_b)
    
    try:
        """
        When accessing the tuples, remember:
            x[0]: index
            x[1]: start
            x[2]: end
        """
        while True:
            """ If a is behind b, spool a """
            
            while a[2] < b[1] or b[2] < a[1]:
                if a[2] < b[1]:
                    a = next(en_a)
                if b[2] < a[1]:
                    b = next(en_b)
            
            region_membership[a[0]] = b[0]
            
            if a[1] < b[1]:
                a = next(en_a)
            else:
                b = next(en_b)
    
    except StopIteration:
        pass
    
    return region_membership


def intersect_ranges_by_chromosome(ranges_a, ranges_b, a_chrom_key="chrom", b_chrom_key="chrom", nonevalue=-1):
    chromosome_intersection = set(ranges_a[a_chrom_key]).intersection(set(ranges_b[b_chrom_key]))
    ranges_a = ranges_a.groupby(a_chrom_key)
    ranges_b = ranges_b.groupby(b_chrom_key)
    
    region_membership = []
    for chrom in chromosome_intersection:
        chr_a = ranges_a.get_group(chrom).drop(a_chrom_key, axis=1)
        chr_b = ranges_b.get_group(chrom).drop(b_chrom_key, axis=1)
        region_membership.append(intersect_ranges(chr_a, chr_b, nonevalue=nonevalue))
    
    return pd.concat(region_membership)
