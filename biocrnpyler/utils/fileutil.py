# findfile.py - find biocrnypyler file
# RMM (via Claude), 19 Sep 2025

import os
from pathlib import Path
from typing import Optional


def find_file_in_bcp_path(
        filename: str, env_var_name: str = "BCP_PATH") -> Optional[str]:
    """
    Find a file by searching through biocrnpyler paths.
    
    Search order:
    1. Current directory
    2. Directories specified in environment variable (default: BCP_PATH)
    3. biocrnpyler package's 'defaults' subdirectory
    
    Args:
        filename (str): Name of file to find
        env_var_name (str): Environment variable containing search paths
        
    Returns:
        str: Full path to file if found, None if not found
        
    Example:
        >>> # Set environment variable (in shell or programmatically)
        >>> os.environ['BCP_PATH'] = '/path/to/models:/path/to/configs'
        >>> 
        >>> # Find a file
        >>> filepath = find_file_in_bcp_path('model.xml')
        >>> if filepath:
        >>>     print(f"Found file at: {filepath}")
        >>> else:
        >>>     print("File not found")
    """
    search_paths = []
    
    # 1. Current directory
    search_paths.append(Path.cwd())
    
    # 2. Environment variable paths
    env_paths = os.environ.get(env_var_name, "")
    if env_paths:
        # Split by colon on Unix, semicolon on Windows
        separator = ";" if os.name == "nt" else ":"
        for path_str in env_paths.split(separator):
            path_str = path_str.strip()
            if path_str:
                path_obj = Path(path_str)
                if path_obj.exists() and path_obj.is_dir():
                    search_paths.append(path_obj)
    
    # 3. biocrnpyler package defaults directory
    import biocrnpyler
    search_paths.append(Path(biocrnpyler.__file__).parent)
    
    # Search for the file in each directory
    for search_dir in search_paths:
        file_path = search_dir / Path(filename)
        if file_path.exists() and file_path.is_file():
            return str(file_path.resolve())

    return None
