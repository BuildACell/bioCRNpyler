# fileutil.py - file utilities for biocrnpyler
# RMM (via Claude), 19 Sep 2025
#
# Utilities for finding and managing files in the biocrnpyler package.

import os
from pathlib import Path
from typing import Optional


def find_file_in_bcp_path(
    filename: str, env_var_name: str = 'BCP_PATH'
) -> Optional[str]:
    """Find a file by searching through biocrnpyler paths.

    Searches for a file in multiple locations in the following order:

    1. Current working directory
    2. Directories specified in environment variable (default: 'BCP_PATH')
    3. biocrnpyler package directory

    The environment variable should contain a colon-separated (Unix) or
    semicolon-separated (Windows) list of directory paths.

    Parameters
    ----------
    filename : str
        Name of file to find.
    env_var_name : str, default='BCP_PATH'
        Name of environment variable containing search paths.

    Returns
    -------
    str or None
        Full path to file if found, None otherwise.

    Examples
    --------
    Set environment variable and find a file:

    >>> import os
    >>> os.environ['BCP_PATH'] = '/path/to/models:/path/to/configs'
    >>> filepath = find_file_in_bcp_path('model.xml')
    >>> if filepath:
    ...     print(f'Found file at: {filepath}')
    ... else:
    ...     print('File not found')
    Found file at: /path/to/models/model.xml

    Search with a custom environment variable:

    >>> filepath = find_file_in_bcp_path(
    ...     'parameters.csv', env_var_name='MY_PARAMS_PATH')

    """
    search_paths = []

    # 1. Current directory
    search_paths.append(Path.cwd())

    # 2. Environment variable paths
    env_paths = os.environ.get(env_var_name, '')
    if env_paths:
        # Split by colon on Unix, semicolon on Windows
        separator = ';' if os.name == 'nt' else ':'
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
