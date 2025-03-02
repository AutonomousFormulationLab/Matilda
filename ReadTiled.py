from tiled.client import from_uri
from tiled.utils import tree
from tiled.queries import FullText


client = from_uri("http://0.0.0.0:8020")

# List available catalogs
catalog_names = list(client)
#print("Available catalogs:", catalog_names)

# Replace 'your_catalog_name' with the actual name of the catalog you want to access
catalog = client[catalog_names[0]]

# List all keys in the catalog
all_keys = list(catalog)
#print("Catalog keys:", all_keys)

# Get the last 10 records
last_10_records = all_keys[-10:]
#print("Last 10 records:", last_10_records)
# Retrieve the last 10 entries using the keys
last_10_entries = [catalog[key] for key in last_10_records]
# Extract titles
#c =  catalog[last_10_records]
for c   in last_10_entries:
    #print(c.metadata['title'])
    titleCont=c.search(FullText("title`"))
    print(titleCont[0])


#LIstOftitles = [entry.metadata['title'] for entry in last_10_entries]
#print("Filenames of the last 1 records:", LIstOftitles)

# List all keys in the catalog
#all_keys = list(catalog)

# Get the last 10 keys
#last_10_keys = all_keys[-10:]

# Retrieve the last 10 entries using the keys
#last_10_entries = [catalog[key] for key in last_10_keys]

#print(last_10_entries)

# Extract UUIDs (assuming each entry has a 'metadata' attribute with a 'uuid' key)
#u#uids = [entry.metadata['uuid'] for entry in last_10_entries]
#print("UUIDs of the last 10 records:", uuids)