import requests
import json
import sys


def nextclade_get_index():
    """Get Nextclade's v3/index.json file and extract the parts useful for matching the
    user's search term and listing available clade annotations.  If the request fails,
    warn and return an empty list."""
    url = "https://data.clades.nextstrain.org/v3/index.json"
    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"Error fetching Nextclade index: {e}", file=sys.stderr)
        return []

    try:
        index_data = response.json()
    except json.JSONDecodeError:
        print("Error decoding JSON from Nextclade index.", file=sys.stderr)
        return []

    entries = []

    for collection in index_data.get('collections', []):
        for dataset in collection.get('datasets', []):
            cladeCount = dataset.get('capabilities', {}).get('clades', 0)
            # Skip entries without clade annotations
            if cladeCount > 0:
                entry = {}
                entry["path"] = dataset.get('path')
                entry["name"] = dataset.get('attributes', {}).get('name')
                entry["segment"] = dataset.get('attributes', {}).get('segment')
                entry["clades"] = {"clade": cladeCount}
                # If customClades are included, break out their counts too
                customClades = dataset.get('capabilities', {}).get('customClades', {})
                if customClades:
                    for customClade, count in customClades.items():
                        entry["clades"][customClade] = count
                entries.append(entry)

    return entries
