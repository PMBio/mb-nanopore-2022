import itertools

from multiprocessing.pool import Pool


class AdvancedPool(Pool):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def __exit__(self, *args, **kwargs):
        super().__exit__(*args, **kwargs)
    
    def __enter__(self, *args, **kwargs):
        super().__enter__(*args, **kwargs)
        return self
    
    def chunky_map_async(self, func, iterable, chunksize, **kwargs):
        """
        Works similar to map_async, but instead of collecting all results in memory and then calling
        the callback once, we call the callback every time a chunk has been completed.
        """
        iterator = iter(iterable)
        jobs = []
        while True:
            chunk = list(itertools.islice(iterator, chunksize))
            if len(chunk) == 0:
                break
            jobs.append(super().map_async(func, chunk, chunksize=chunksize, **kwargs))
        for job in jobs:
            job.wait()
    
    def chunky_starmap_async(self, func, iterable, chunksize, **kwargs):
        """
        Works similar to starmap_async, but instead of collecting all results in memory and then calling
        the callback once, we call the callback every time a chunk has been completed.
        """
        iterator = iter(iterable)
        jobs = []
        while True:
            chunk = list(itertools.islice(iterator, chunksize))
            if len(chunk) == 0:
                break
            jobs.append(super().starmap_async(func, chunk, chunksize=chunksize, **kwargs))
        for job in jobs:
            job.wait()
