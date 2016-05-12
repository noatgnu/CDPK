import urllib.request
import urllib.parse
import csv
import re
import io
from uniprot_parser import settingup
import os
from tornado import gen
from tornado.httpclient import AsyncHTTPClient
from tornado.ioloop import IOLoop
from redis import Redis
from rq import Queue
import time
def allowed_file(filename):
    split_fn = filename.split('.')[-1]
    if split_fn:
       if split_fn in settingup.ALLOWED_EXTENSIONS:
           return 'True'
       else:
           return 'False'
    return 'False'

def check_column(input_path):
    accession_col_name = ''
    with open(input_path, 'rt', newline='') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel', delimiter='\t')
        for row in reader:
            for k in row:
                if re.search('\|[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\|', str(row[k])):
                    accession_col_name = k
            break
    return accession_col_name

def file_set_up(input):
    input_path = os.path.join(settingup.UPLOAD_FOLDER, input)
    filen = input.replace('.txt', '').replace('.csv', '')
    return input_path, filen


def uniprot_parser(protid, infocolumn):
    user_agent = "toan.phung@uqconnect.uq.edu.au"
    headers = {'User Agent': user_agent}
    contact = {'name': 'Toan', 'Institution': 'University of Queensland', 'Language': 'Python'}

    base_url = 'http://www.uniprot.org/uniprot/?query=accession:'
    query = base_url+protid+"&"+"format=tab"+'&'+'compression=no'+'&'+'columns='+infocolumn
    
    stuff = urllib.parse.urlencode(contact)
    stuff = stuff.encode('utf-8')
    #http_client = AsyncHTTPClient()
    req = urllib.request.Request(query, stuff, headers)
    response = urllib.request.urlopen(req)
    #response = yield http_client.fetch(query)
    response = response.read()
    response = response.decode('utf-8')
    data = io.StringIO(response)
    csv_data = csv.DictReader(data, delimiter='\t')
    parsing_data = {}
    
    for row in csv_data:
        for k in row:
            if k == 'Subcellular location [CC]':
                v = row[k].replace('SUBCELLULAR LOCATION: ', '')
                parsing_data[k] = v
            else:
                parsing_data[k] = row[k]
        break        
    return parsing_data


def file_processing(input_path, filen, accession_col_name):
    n = 0
    OutputFile = os.path.join(settingup.RESULTS_FOLDER, filen+'matched_uniprot.txt')
    OutputFile2 = os.path.join(settingup.RESULTS_FOLDER, filen+'matched_uniprot_gs.txt')
    OutputFile22 = os.path.join(settingup.RESULTS_FOLDER, filen+'matched_uniprot_goen.txt')
    OutputFile3 = os.path.join(settingup.RESULTS_FOLDER, filen+'matched_uniprot_gs_no_goterm.txt')
    
    field_names = ['Entry', 'Entry name', 'Protein names', 'Subcellular location [CC]', 'Gene ontology (GO)', 'Gene ontology IDs', 'Sequence', 'N-glycosylation Sequon Positions', 'Glycosylation']
    field_names_gs = ['Goterm', 'Gocode', 'Entry']
    field_names_goen = ['Entry', 'GoIDs']
    with open(input_path, 'rt', newline='') as csvfile, open(OutputFile, 'wt', newline='', encoding='utf-8') as resource, open(OutputFile2, 'wt', newline='', encoding='utf-8') as resource2, open(OutputFile22, 'wt', newline='', encoding='utf-8') as resource22,open(OutputFile3, 'wt', newline='', encoding='utf-8') as resource3:
        reader = csv.DictReader(csvfile, dialect='excel', delimiter='\t')
        output_write = csv.DictWriter(resource, fieldnames=field_names, dialect='excel', delimiter='\t')
        output_write.writeheader()
        output_write2 = csv.DictWriter(resource2, fieldnames=field_names_gs, dialect='excel', delimiter='\t')
        output_write2.writeheader()
        output_write22 = csv.DictWriter(resource22, fieldnames=field_names_goen, dialect='excel', delimiter='\t')
        output_write22.writeheader()
        output_write3 = csv.DictWriter(resource3, fieldnames=field_names_gs, dialect='excel', delimiter='\t')
        output_write3.writeheader()

        n = 0
        for row in reader:
            n += 1
            
            gsoutput = {}
            goenoutput = {}
            rfilter = re.search('\|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\|', row[accession_col_name])
            #query_str = query_str + '+' + rfilter.group(1)
            #if n == 3:
               # n = 0                
                #print(query_str)
            uniprot_data = uniprot_parser(rfilter.group(1), 'id,entry%20name,protein%20names,comment(SUBCELLULAR%20LOCATION),go,go-id,sequence,feature(GLYCOSYLATION)')
            uniprot_data['N-glycosylation Sequon Positions'] = ''
            if not uniprot_data['Sequence'] == '':
                pos_search = re.finditer("N\w(S|T)", uniprot_data['Sequence'])
                position = ''
                for i in pos_search:
                    if not i.group(0)[1] == 'P':
                        position = position + str(i.start()) + '-' + str(i.end()) + ';'
                uniprot_data['N-glycosylation Sequon Positions'] = position
            gsoutput['Entry'] = uniprot_data['Entry name']
            gsoutput['Gocode'] = 'IEA'
            if uniprot_data['Gene ontology IDs'] == '':
                gsoutput['Goterm'] = ''
                output_write3.writerow(gsoutput)
            else:
                goenoutput['Entry'] = uniprot_data['Entry name']
                goterm_ids_string = uniprot_data['Gene ontology IDs'].replace(" ", "")
                goenoutput['GoIDs'] = goterm_ids_string
                output_write22.writerow(goenoutput)
                goterm_id_list = goterm_ids_string.split(';')
                for i in goterm_id_list:
                    gsoutput['Goterm'] = i
                    output_write2.writerow(gsoutput)
            output_write.writerow(uniprot_data)
                                                                                                                          
    return {filen: [filen+'matched_uniprot.txt', filen+'matched_uniprot_gs.txt', filen+'matched_uniprot_goen.txt', filen+'matched_uniprot_gs_no_goterm.txt']}


def processor(filen): 
    input_fp, fn = file_set_up(filen)
    acc_col_name = check_column(input_fp)
    output = file_processing(input_fp, fn, acc_col_name)
    return output


def get_results(job_key):
    q = Queue(connection=Redis())
    j = q.fetch_job(job_key)
    job_status = ''
    while True:
        time.sleep(30) 
        if j.status == 'finished':
            job_status = 'finished'
            break
        if j.status == 'failed':
            job_status = 'failed'
            break
    if job_status == 'finished':         
        return j.result
    if job_status == 'failed':
        return {'Failed': ['Job failure. Please make sure that the file format was correct.']}
