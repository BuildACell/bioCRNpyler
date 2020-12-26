#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from pathlib import Path
import json

IGNORE = [
    "3. DNA Assemblies gene expression transcription and translation",
     "4. Promoters Transcriptional Regulation and Gene Regulatory Networks",
]
# offset for the nblink files inside source/examples
offset = '../../'


def create_nblink(notebook_path, offset):
    d = {}
    d["path"] = f'{offset}{notebook_path}'
    # d["extra-media"] = [str(notebook_path.parent / "images/")]
    return d


def main():
    examples = Path("../examples/")
    list_of_examples = []
    files = sorted(examples.glob("*.ipynb"))
    for one_nb in files:

        if not one_nb.exists():
            continue
        for ignore in IGNORE:
            if ignore in str(one_nb):
                break
        else:
            nblink = create_nblink(one_nb, offset)
            list_of_examples.append(f'.. toctree::\n \t:maxdepth: 2\n \t:caption: {one_nb.stem}\n\n \texamples/{one_nb.stem}.ipynb\n\n')
            with open(f"source/examples/{one_nb.stem}.nblink", "w") as f:
                print("Creating path for", one_nb)
                json.dump(nblink, f)

    # write nb_examples.rst
    with open(f"source/nb_examples.rst", "w") as f:
        print("Creating nb_examples.rst")
        [f.writelines(l) for l in list_of_examples]


if __name__ == "__main__":
    main()