# readTiled.py 
# version 0.1   2025-03-04


# this code will get last N data sets from tiled database
# it will recreate what IR3BS_GetJSONScanData() does in Indra package In3_ClueSkyreader.ipf
# it wil luse tiled web interface and read data into json.
# the main purpose is to get list fo filenames and paths for the last N scans
# and then use this list to download the data from the server
#the code will aslo remeber when was the last polling and ask for data only since last time. 
# need to handel gracefully the start
# code will trigger appropriate data recduction routines for different data sets
# at the end the code will generate necessary pictures for live data page.
# the code will be run in background and be checked by a cron job every 5 minutes

# method used buitls on https://github.com/BCDA-APS/bdp-tiled/blob/main/demo_client.ipynb
# and follows Igor code to get the right data sets


# import necessary libraries
import requests
import json
import datetime


def iso_to_ts(isotime):
    return datetime.datetime.fromisoformat(isotime).timestamp()

server = "localhost"
port = 8020
catalog = "usaxs"

# Find all runs in a catalog between these two ISO8601 dates.
start_time = "2025-02-15"
end_time = "2025-02-25"
tz = "US/Central"


#  http://10.211.55.7:8020/api/v1/search//usaxs/?page[offset]=0
# &page[limit]=100
# &filter[time_range][condition][since]=1738738800
# &filter[time_range][condition][timezone]=US/Central
# &sort=time
# &fields=metadata
# &omit_links=true
# &select_metadata={detectors:start.detectors,motors:start.motors,plan_name:start.plan_name,time:start.time,
# scan_title:start.plan_args.scan_title,title:start.title,hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}
uri =(
    f"http://{server}:{port}"  # standard prefix
    "/api/v1/search"    # API command
    f"/{catalog}"       # catalog
    "?"                 # begin any command options
    "page[limit]=10"   # 0: all matching
    "&"                 # separator between any additional options
    f"filter[time_range][condition][since]={iso_to_ts(start_time)}"
    f"&filter[time_range][condition][until]={iso_to_ts(end_time)}"
    f"&filter[time_range][condition][timezone]={tz}"
    "&sort=time"
    "&fields=metadata"
    "&omit_links=true"
    "&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
                        hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"
)
#print(f"{uri=}")
r = requests.get(uri).json()
#print(r["data"][0])

def print_results_summary(r):
    """We'll use this a few times."""
    xref = dict(First=0, Last=-1)
    for k, v in dict(First=0, Last=-1).items():
        md = r["data"][v]["attributes"]["metadata"]["selected"]  #this is for UBUNTU VM, usaxscontrol does not have selected
        #print(md)
        # md keys: start  stop  summary
        # summary key is composed by tiled server
        plan_name = md["plan_name"]
        #scan_id = md["data"]["id"]
        #started = md["data"]["datetime"]
        hdf5_file = md["hdf5_file"]
        hdf5_path = md["hdf5_path"]
        print(f"{k:5s} run: {plan_name=} path: {hdf5_path=} {hdf5_file=}")
        
print(f'Search of {catalog=} has {len(r["data"])} runs.')
print_results_summary(r)


