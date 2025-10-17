import csv
import json
import os
import time
from functools import partial
import boto3
from tqdm.contrib.concurrent import thread_map

def restore_object(key, bucket, tier):
    """
    Check Object status before requesting the restoration.
    returns a dict with key and the status of the request. If there is an error
    message and object metadata is also included in the response.
    """
    client = boto3.session.Session().client("s3")
    metadata = client.head_object(Bucket=bucket, Key=key)
    headers = metadata["ResponseMetadata"]["HTTPHeaders"]
    if "x-amz-archive-status" not in headers:
        return {"key": key, "status": "RESTORED"}
    if "Restore" in metadata:
        restore_tag = metadata["Restore"]
        if "expiry-date" in restore_tag:
            return {"key": key, "status": "RESTORED"}
        if 'ongoing-request="true"' in restore_tag:
            return {"key": key, "status": "IN_PROGRESS"}
    try:
        client.restore_object(Bucket=bucket, Key=key, RestoreRequest={'GlacierJobParameters':{'Tier':tier}})
        return {"key": key, "status": "REQUESTED"}
    except Exception as ex:
        return {"key": key,
                "status": "ERROR",
                "message": {str(ex)},
                "metadata": json.dumps(metadata, default=str)}

def summarize_and_log_results(outputs,logfile):
    """
    Summarize and log to a file the restoration results produced
    by running restore_object. Return the errors, in case we want to
    try again.
    """
    error_list = []
    counts = {'REQUESTED': 0, 'IN_PROGRESS': 0, 'RESTORED': 0, 'ERROR': 0}
    if os.path.dirname(logfile)!= '':
        if not os.path.exists(os.path.dirname(logfile)):
            os.makedirs(os.path.dirname(logfile),exist_ok=True)
    with open(logfile, 'w', newline='') as csvfile:
        writer = csv.DictWriter(
                csvfile, fieldnames=["key", "status", "message", "metadata"])
        writer.writeheader()
        for output in outputs:
            counts[output["status"]] += 1
            writer.writerow(output)
            if output["status"] == "ERROR":
                error_list.append(output["key"])
    for status, count in counts.items():
        print(f'{status:<14} {count}')
    print(f'For more info check {logfile}')
    return error_list

def bulk_restore(bucket,prefix,is_logfile=False,filter_in=None,filter_out=None,tier='Standard',max_workers=8,logfile='output.csv',retry_once=False):
    """
    Bulk restore a bunch of things in S3 that are in the IntelligentTiering class.
    You need to pass in a bucket and either a) a prefix (aka file OR folder) or b) a file name (see below).
    If your prefix is a folder, you can optionally then also pass some filters to subset the data.
    filter_in and filter_out can both take either strings or lists.
    If a list in filter_in, it will keep any object that has ANY filter.
    If a list in filter_out, it will remove any object that has ANY filter.
    You should always use tier=Standard for large requests unless it is a very severe emergency,
    (A standard cell painting plate is ~$300 to restore in Expedited but free in Standard)
    but if it IS an emergency or it's only a couple files, you can pass tier='Expedited' instead).

    Sometimes, errors will happen!
    You can optionally use --retry_once to go through the list of failed restorations at the end of your run.
    This is good for small one off "hiccups", but doing it right away may be too soon though to try again
    if large numbers of things failed (ie if you're throttled). You can therefore also optionally pass in
    the name of the log file from a previous run in order to try the errors again. Passing in a file name
    will override any prefix set.

    To call at the command line, use `python restore_intelligent.py bucketname folder --filter_in ch1 ch2 --filter_out .txt .csv`
    """
    file_list = []
    if not is_logfile:
        client = boto3.client('s3')
        other_tier_count = 0
        paginator = client.get_paginator("list_objects_v2")
        pages = paginator.paginate(Bucket=bucket, Prefix=prefix)
        try:
            for page in pages:
                file_list += [x["Key"] for x in page["Contents"] if x["StorageClass"] == "INTELLIGENT_TIERING"]
                other_tier_count += len([x["Key"] for x in page["Contents"] if x["StorageClass"] != "INTELLIGENT_TIERING"])
        except KeyError:
            print ("No files in prefix given.")
            return
        print(f"{other_tier_count} non-Intelligent Tiering files that do not need restoration")
        print(f"{len(file_list)} total files found pre-filtering")
    else:
        count = -1
        if os.path.exists(prefix):
            with open(prefix) as csvfile:
                reader = csv.reader(csvfile)
                for row in reader:
                    count+=1
                    if row[1] == 'ERROR':
                        file_list.append(row[0])
            print(f"{len(file_list)} previous errors found in {count} rows checked in {prefix}")
        else:
            print(f"{prefix} not found")
            return
    if filter_in:
        if type(filter_in) == str:
            file_list = [x for x in file_list if filter_in in x]
        elif type(filter_in) == list:
            masterlist = []
            for eachfilter in filter_in:
                filtered_list = [x for x in file_list if eachfilter in x]
                masterlist += filtered_list
            file_list = list(set(masterlist))
    if filter_out:
        if type(filter_out) == str:
            file_list = [x for x in file_list if filter_out not in x]
        elif type(filter_out) == list:
            for eachfilter in filter_out:
                file_list = [x for x in file_list if eachfilter not in x]
    print(f"{len(file_list)} total files remain post-filtering")
    file_list.sort() #It is nice to have this in order in case you need to try to manually figure out where it left off/how far it has gotten
    outputs = thread_map(partial(restore_object, tier=tier, bucket=bucket),
                         file_list, max_workers=max_workers)
    error_list = summarize_and_log_results(outputs,logfile)
    if retry_once and len(error_list) > 0:
        print(f"Retrying {len(error_list)} errored files")
        retry_outputs = thread_map(partial(restore_object, tier=tier, bucket=bucket),
                     error_list, max_workers=max_workers)
        summarize_and_log_results(retry_outputs,logfile[:-4]+'_retried.csv')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Restore Intelligent-Tiering objects in bulk')
    parser.add_argument('bucket', help='Name of the bucket')
    parser.add_argument('prefix', help='Prefix (aka file or folder) to restore. If --is_logfile is passed, this should instead be a logfile from a previous run')
    parser.add_argument('--is_logfile',action="store_true",default=False,help='Indicates the prefix is a log file from a previous run; instead listing objects in the bucket, this will be searched for errored keys and those will be retried')
    parser.add_argument('--filter_in',default=None,nargs='+',help='One or more strings to specify a subset of objects to restore. If >1, any object with ANY filter will be restored')
    parser.add_argument('--filter_out',default=None,nargs='+',help='One or more strings to specify a subset of objects not to restore. If >1, any object with ANY filter will be removed from the list to restore')
    parser.add_argument('--tier',default='Standard',help='Retrieval tier, only change to Expedited for small numbers of files (<1K) or emergencies')
    parser.add_argument('--max_workers',default=8,help='Number of parallel AWS requests', type=int)
    parser.add_argument('--logfile',default='output.csv',help='Path to save the status log in csv format')
    parser.add_argument('--retry_once',action="store_true",default=False,help='Optionally retried failed files one more time')
    args = parser.parse_args()

    bulk_restore(args.bucket,args.prefix,args.is_logfile,args.filter_in,args.filter_out,
                 args.tier,args.max_workers,args.logfile,args.retry_once)
