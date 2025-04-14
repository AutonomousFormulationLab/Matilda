'''

notes on tiled communication

'''


def FindAscan(plan_name,matchStr, NumScans=1):
    #this filters for last collected scan matching specific sttring in title 
    uri = (
        f"http://{server}:{port}"
        "/api/v1/search"
        f"/{catalog}"
        f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
        "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
        f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
        "&filter[regex][condition][key]=title"                              # filter by title
        f'&filter[regex][condition][pattern]=(?i){matchStr}'                   # filter - using gtrep - for title matchStr, casxe independent
        "&sort=-time"                                                       # sort by time, -time gives last scans first
        "&fields=metadata"                                                  # return metadata
        "&omit_links=true"                                                  # no links
        "&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
        hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"   # select metadata
        )
    logging.info(f"{uri=}")
    #print(f"{uri=}")
    #additional methods to match data:
    #&filter[regex][condition][key]=plan_name
    #&filter[regex][condition][pattern]={plan_name}     #use regex to match plan_name
    #filter[regex][condition][key]=title
    #&filter[regex][condition][key]={title}             #use regex to match title  
    #and
    #[noteq] - not equal
    #[contains] - seems same as eq in use, and cannot be made into case insensitive. Not useful. 
    #[in] - in a list of values
    #[notin] - not in a list of values
    #[comparison] - comparison with lt, gt, le, ge for numerical values
    #working example:
    #http://10.211.55.7:8020/api/v1/search/usaxs/?page[limit]=10&filter[eq][condition][key]=plan_name&filter[eq][condition][value]=%22WAXS%22&filter[regex][condition][key]=title&filter[regex][condition][pattern]=(?i)blank&sort=-time
    #returns list of "Blank" samples, not not ist of samples contains "blank" in name
    #http://10.211.55.7:8020/api/v1/search/usaxs/?page[limit]=1&filter[eq][condition][key]=plan_name&filter[eq][condition][value]=%22WAXS%22&filter[regex][condition][key]=title&filter[regex][condition][pattern]=(?i)water*blank&sort=-time
    #returns last scan which conatins case independent "water blank" in name
    #http://10.211.55.7:8020/api/v1/search/usaxs/?page[limit]=1&filter[eq][condition][key]=plan_name&filter[eq][condition][value]=%22WAXS%22&filter[regex][condition][key]=title&filter[regex][condition][pattern]=(?i)blank&sort=-time&omit_links=true&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}
    #returns last scan which conatisn case independet "water blank" in name