import shutil
import os
from pathlib import Path
from typing import List

def get_top_existing_dir(path: Path) -> Path:
    if path.exists():
        return path
    else:
        return get_top_existing_dir(path.parent)

def get_free_space_gb(path: str) -> float:
    path = Path(path).resolve()
    path = get_top_existing_dir(path)
    gb_free = shutil.disk_usage(path).free / (1024 * 1024 * 1024)
    return gb_free

def has_enough_free_space(path, gb_required: float) -> bool:
    gb_free = get_free_space_gb(path)
    print(f"Need {gb_required} has {gb_free}")
    return gb_free > gb_required

def check_free_space(path: str, gb_required: float) -> str:
    if has_enough_free_space(path, gb_required):
        return path
    raise RuntimeError(f"Not enough disk space on {path}. Need {gb_required}gb.")

def compute_total_file_size_gb(file_list: List[str]) -> float:
    return sum(Path(f).stat().st_size for f in file_list) / (1024*1024*1024)
