# Description: This code snippet demonstrates how to find the path to a specific key in a JSON file.
# The find_key_path function recursively searches through the JSON data to find the specified key.
# If the key is found, the function returns the path to the key as a list of keys and indices.
# If the key is not found, the function returns None.
# The script loads a JSON file, specifies the key to find, and calls the find_key_path function to search for the key.

import json

def find_key_path(data, target_key, path=None):
    if path is None:
        path = []

    if isinstance(data, dict):
        for key, value in data.items():
            new_path = path + [key]
            if key == target_key:
                return new_path
            elif isinstance(value, (dict, list)):
                result = find_key_path(value, target_key, new_path)
                if result is not None:
                    return result

    elif isinstance(data, list):
        for index, item in enumerate(data):
            new_path = path + [index]
            result = find_key_path(item, target_key, new_path)
            if result is not None:
                return result

    return None

# Load your JSON file
with open('your_file.json', 'r') as file:
    json_data = json.load(file)

# Specify the key you are looking for
key_to_find = 'your_target_key'

# Find the path to the key
path = find_key_path(json_data, key_to_find)

if path is not None:
    print(f"Path to '{key_to_find}': {path}")
else:
    print(f"Key '{key_to_find}' not found in the JSON data.")
