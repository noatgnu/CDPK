import tornado.ioloop
import tornado.web
from sequon import settingsequon, seqtasks
from uniprot_parser import settingup, uptasks, uniparser
import os
from tornado.httpclient import AsyncHTTPClient
from tornado import gen, websocket, httpserver
import settingmain
from concurrent.futures import ThreadPoolExecutor
from redis import Redis
from rq import Queue
import json
import time
from rq.job import Job
from tornado.escape import json_encode
from csvcompare import comparesetting
from goenrichment import gosetting
import csv
import MySQLdb
import re
import subprocess
import rpy2.robjects as robjects 
from goenrichment import commonf

class wsUI(tornado.web.RequestHandler):
    def get(self):
        self.render("websocket.html")
        
class WebSocketHandler(websocket.WebSocketHandler):
    def open(self):
        pass
        # Websocket function for Uniprot Parser version that requires an internet connection.
    @gen.coroutine
    def on_message(self, parameters):
        input_p = json.loads(parameters) # Decode JSON dictionary into a python dictionary to get key information.
        pool = ThreadPoolExecutor(max_workers=settingmain.MAX_WORKERS)
        if input_p['task'] == 'uniprotparser':
            # Input key to look for result.
            result = yield pool.submit(uptasks.get_results, input_p['key'])
            pool.shutdown()
            self.write_message(result)
        if input_p['task'] == 'uniprotparser2':
            # This statement have no function as of yet.
            result = yield pool.submit(uniparser.uniprot_parser, input_p)
            pool.shutdown()
            self.write_message(json.encode(result))
         
    def on_close(self):
        pass


class MainHandler(tornado.web.RequestHandler):
    def get(self):
        self.render("index.html", title="Home Page",)
        # Index page.
        
class SequonParserUploadForm(tornado.web.RequestHandler):
    def get(self):
        self.render("upload_sequon.html", title="Upload Peptide Summary",)
        # Upload page for sequon parser.
        
class SequonParserUploadFileProcess(tornado.web.RequestHandler):
    @gen.coroutine
    def post(self):
        # Processing the uploaded file
        fileinfo = self.request.files['file'][0]
        filen = fileinfo['filename']
        if seqtasks.allowed_file(filen) == 'True':
            # If the file's extension is .txt or .csv, store the file.
            with open(os.path.join(settingsequon.UPLOAD_FOLDER, filen), 'wb') as upload_file:
                upload_file.write(fileinfo['body'])
            # Execute the function on the just stored file.
            pool = ThreadPoolExecutor(max_workers=settingmain.MAX_WORKERS)
            result_fn, excel = yield pool.submit(seqtasks.processor, filen)
            pool.shutdown()
            # Render result.
            self.render("results_sequon.html", title="Results", results = result_fn, excel_list = excel,)

class UniprotUploadForm(tornado.web.RequestHandler):
    def get(self):
        self.render("upload_uniprot.html", title="Upload SeqID",)
        # Upload page for Uniprot Parser, Online database version.
        
class UniprotUploadFileProcess(tornado.web.RequestHandler):
    @gen.coroutine
    def post(self):
        fileinfo = self.request.files['file'][0]
        filen = fileinfo['filename']
        if uptasks.allowed_file(filen) == 'True':
            # Check for allowed extensions and store the file.
            with open(os.path.join(settingup.UPLOAD_FOLDER, filen), 'wb') as upload_file:
                upload_file.write(fileinfo['body'])
            
            # Move the job to rq and off main thread.
            q = Queue('uniprotparser', connection=Redis())
            # Job will timeout after 14400 seconds if has not finished yet.
            result = q.enqueue(uptasks.processor, args=(filen,), timeout=14400)
            # Dictionary for storing result
            output = {'task':'uniprotparser', 'key':result.id}
            # Render result.
            self.render("results_uniprot.html", title="Results", result=output,)
            
class GetResult(tornado.web.RequestHandler):
    def get(self, job_key, channel_id):
        
        q = Queue(channel_id, connection=Redis())
        j = q.fetch_job(job_key)
        if j.status == 'finished':
            self.write(json_encode(j.result))
        elif j.status == 'failed':
            self.write('Fail')
        else:
            self.set_status(202)
            self.write("Unfinished")
class CSVCompareForm(tornado.web.RequestHandler):
    def get(self):
        #CSVComparison Upload site
        self.render("compare.html", title="Compare",)
        
class CSVCompareProcess(tornado.web.RequestHandler):
    
    def post(self):
        fileinfo1 = self.request.files['file1'][0]
        fileinfo2 = self.request.files['file2'][0]
        filen1 = fileinfo1['filename']
        filen2 = fileinfo2['filename']
        
        # Check uploaded files' extension and save files, protein summary only.
        if seqtasks.allowed_file(filen1) == 'True' and seqtasks.allowed_file(filen2) == 'True':
            with open(os.path.join(comparesetting.UPLOAD_FOLDER, filen1), 'wb') as upload_file1, open(os.path.join(comparesetting.UPLOAD_FOLDER, filen2), 'wb') as upload_file2:
                upload_file1.write(fileinfo1['body'])
                upload_file2.write(fileinfo2['body'])
            acc_list1 = {}
            acc_list2 = {}
            # Generate dictionaries with keys from the accession ids of each file.
            with open(os.path.join(comparesetting.UPLOAD_FOLDER, filen1), 'rt', newline='') as csvfile1, open(os.path.join(comparesetting.UPLOAD_FOLDER, filen2), 'rt', newline='') as csvfile2: 
                reader1 = csv.DictReader(csvfile1, dialect='excel', delimiter='\t')
                reader2 = csv.DictReader(csvfile2, dialect='excel', delimiter='\t')
                n = 0
                m = 0
                for row in reader1:
                    n += 1
                    acc_list1[row['Accession']] = n
                for row in reader2:
                    m += 1
                    acc_list2[row['Accession']] = m
            # Detect differences in the two dictionary and write them out in different files.        
            with open(os.path.join(comparesetting.UPLOAD_FOLDER, filen1), 'rt', newline='') as csvfile1, open(os.path.join(comparesetting.UPLOAD_FOLDER, filen2), 'rt', newline='') as csvfile2, open(os.path.join(comparesetting.RESULTS_FOLDER, 'not_in_'+filen1), 'wt', newline='', encoding='utf-8') as output1, open(os.path.join(comparesetting.RESULTS_FOLDER, 'not_in_'+filen2), 'wt', newline='', encoding='utf-8') as output2:             
                reader1 = csv.DictReader(csvfile1, dialect='excel', delimiter='\t')
                reader2 = csv.DictReader(csvfile2, dialect='excel', delimiter='\t')
                field_names2 = reader2.fieldnames
                field_names2 = list(filter(None, field_names2))
                field_names1 = reader1.fieldnames
                field_names1 = list(filter(None, field_names1))
                #print(field_names2)
                output1_write = csv.DictWriter(output1, fieldnames=field_names2, dialect='excel')
                output1_write.writeheader()
                output2_write = csv.DictWriter(output2, fieldnames=field_names1, dialect='excel')
                output2_write.writeheader()
                for row2 in reader2:
                    row3 = {}
                    if not row2['Accession'] in acc_list1:
                        for i in field_names2:
                            row3[i] = row2[i]
                        output1_write.writerow(row3)
                for row1 in reader1:
                    row3 = {}
                    if not row1['Accession'] in acc_list2:
                        for i in field_names1:
                            row3[i] = row1[i]
                        output2_write.writerow(row3)
            self.render("results_compare.html", title="Compare Result", results=['not_in_'+filen1, 'not_in_'+filen2])

class SubcellLocAnalyze(tornado.web.RequestHandler):
    def get(self):
        self.render("subcell.html", title="Subcellular",)
        
class SubcellLocAnalyzeProcess(tornado.web.RequestHandler):
    
    def post(self):
        proteininfo = self.request.files['file1'][0]
        peptide = self.request.files['file2'][0]
        proteininfo_n = proteininfo['filename']
        peptide_n = peptide['filename']
        if seqtasks.allowed_file(proteininfo_n) == 'True' and seqtasks.allowed_file(peptide_n) == 'True':
            with open(os.path.join(comparesetting.UPLOAD_FOLDER, proteininfo_n), 'wb') as upload_file1, open(os.path.join(comparesetting.UPLOAD_FOLDER, peptide_n), 'wb') as upload_file2:
                upload_file1.write(proteininfo['body'])
                upload_file2.write(peptide['body'])
            acc_er = {}
            acc_golgi = {}
            acc_vacuole = {}
            acc_membrane = {}
            with open(os.path.join(comparesetting.UPLOAD_FOLDER, proteininfo_n), 'rt', newline='') as csvfile1: 
                reader1 = csv.DictReader(csvfile1, dialect='excel', delimiter='\t')                
                for row in reader1:
                    # Build dictionary of proteins and their subcellular locations.
                    # Sorting by order, prioritizing ER, then golgi, then vacuole. If those key words were not found and Mitochondrion was found, ignore entry.
                    # After all conditions above and the entries still not processed and key word 'membrane' was found, write it to new dictionary.
                    if row['Subcellular location [CC]'].find('Endoplasmic reticulum') >= 0:
                        acc_er[row['Entry']] = row['Subcellular location [CC]']
                        continue
                    if row['Subcellular location [CC]'].find('Golgi') >= 0:
                        acc_golgi[row['Entry']] = row['Subcellular location [CC]']
                        continue
                    if row['Subcellular location [CC]'].find('Vacuole') >= 0:
                        acc_vacuole[row['Entry']] = row['Subcellular location [CC]']
                        continue
                    if row['Subcellular location [CC]'].find('Mitochondrion') >= 0:
                        continue
                    if row['Subcellular location [CC]'].find('membrane') >= 0:
                        acc_membrane[row['Entry']] = row['Subcellular location [CC]']
                
            with open(os.path.join(comparesetting.UPLOAD_FOLDER, peptide_n), 'rt', newline='') as csvfile2, open(os.path.join(comparesetting.RESULTS_FOLDER, 'ER_'+peptide_n), 'wt', newline='', encoding='utf-8') as output1, open(os.path.join(comparesetting.RESULTS_FOLDER, 'Golgi_'+peptide_n), 'wt', newline='', encoding='utf-8') as output2, open(os.path.join(comparesetting.RESULTS_FOLDER, 'Vacuo_'+peptide_n), 'wt', newline='', encoding='utf-8') as output3, open(os.path.join(comparesetting.RESULTS_FOLDER, 'Membrane_'+peptide_n), 'wt', newline='', encoding='utf-8') as output4:            # From the accession IDs, match them with the subcellular location dictionaries and write out result accordingly.         
                reader2 = csv.DictReader(csvfile2, dialect='excel', delimiter='\t')
                field_names2 = reader2.fieldnames
                field_names2 = list(filter(None, field_names2))
                field_names2.append('Subcellular location')
                
                output1_write = csv.DictWriter(output1, fieldnames=field_names2, dialect='excel', delimiter='\t')
                output1_write.writeheader()
                output2_write = csv.DictWriter(output2, fieldnames=field_names2, dialect='excel', delimiter='\t')
                output2_write.writeheader()
                output3_write = csv.DictWriter(output3, fieldnames=field_names2, dialect='excel', delimiter='\t')
                output3_write.writeheader()
                output4_write = csv.DictWriter(output4, fieldnames=field_names2, dialect='excel', delimiter='\t')
                output4_write.writeheader()
                for row2 in reader2:
                    
                    rfilter = re.search('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', row2['Accessions'])
                    if rfilter.group(0) in acc_er:
                        row2['Subcellular location'] = acc_er[rfilter.group(0)]
                        output1_write.writerow(row2)
                    if rfilter.group(0) in acc_golgi:
                        row2['Subcellular location'] = acc_golgi[rfilter.group(0)]
                        output2_write.writerow(row2)
                    if rfilter.group(0) in acc_vacuole:
                        row2['Subcellular location'] = acc_vacuole[rfilter.group(0)]
                        output3_write.writerow(row2)
                    if rfilter.group(0) in acc_membrane:
                        row2['Subcellular location'] = acc_membrane[rfilter.group(0)]
                        output4_write.writerow(row2)
                        
            self.render("results_subcell.html", title="Subcellular Analysis Result", results={peptide_n: ['ER_'+peptide_n, 'Golgi_'+peptide_n, 'Vacuo_'+peptide_n, 'Membrane_'+peptide_n]},)
            
class GOEnrichmentForm(tornado.web.RequestHandler):
    def get(self):
        self.render("go.html", title="GO Enrichment",)
        # GO Enrichment file submission.
        
class GOEnrichmentResult(tornado.web.RequestHandler):
    @gen.coroutine
    def post(self):
        parameters = dict()
        # Create dictionary containing all the input parameters.
        studyfile = self.request.files['study'][0] # Study sample
        populationfile = self.request.files['population'][0] # Universe sample
        associationfile = self.request.files['association'][0] # Association
        
        pvd = self.request.body_arguments["pvalue"][0]
        
        parameters['pvalue'] = float(pvd.decode('utf-8'))
        od = self.request.body_arguments["organism"][0]
        
        parameters['organism'] = od.decode('utf-8')
        parameters['ontology'] = []
        if 'mf' in self.request.body_arguments:
            parameters['ontology'].append('MF')
        if 'bp' in self.request.body_arguments:
            parameters['ontology'].append('BP')
        if 'cc' in self.request.body_arguments:
            parameters['ontology'].append('CC')
        parameters['study'] = studyfile['filename']
        parameters['universe'] = populationfile['filename']
        parameters['association'] = associationfile['filename']
        ds = self.request.body_arguments["dstudy"][0]
        parameters['dstudy'] = ds.decode('utf-8')
        du = self.request.body_arguments["duniverse"][0]
        parameters['duniverse'] = du.decode('utf-8')
        da = self.request.body_arguments["dassociation"][0]
        parameters['dassociation'] = da.decode('utf-8')
        
        # Check uploaded files' extension. If correct, save them.
        if seqtasks.allowed_file(parameters['study']) == 'True' and seqtasks.allowed_file(parameters['universe']) == 'True' and seqtasks.allowed_file(parameters['association']) == 'True':
            with open(os.path.join(gosetting.UPLOAD_FOLDER, parameters['study']), 'wb') as upload_file1, open(os.path.join(gosetting.UPLOAD_FOLDER, parameters['universe']), 'wb') as upload_file2, open(os.path.join(gosetting.UPLOAD_FOLDER, parameters['association']), 'wb') as upload_file3:
                upload_file1.write(studyfile['body'])
                upload_file2.write(populationfile['body'])
                upload_file3.write(associationfile['body'])
            pool = ThreadPoolExecutor(max_workers=settingmain.MAX_WORKERS)
            # Transfer user's input to the go enrichment function.
            result = yield pool.submit(commonf.GOEnri, parameters)
            pool.shutdown()
            self.render("results_go.html", title="GO Enrichment Result", results=result,)

class UniparserForm(tornado.web.RequestHandler):
    def get(self):
        self.render("uniparser2.html", title="Uniprot Parser", dboptions=settingup.TABLE_NAME,)
        # Upload form for Uniprotparser-EX
        
class UniparserProcess(tornado.web.RequestHandler):
    @gen.coroutine
    def post(self):
    
        fileinfo = self.request.files['inputfile'][0]
        filen = fileinfo['filename']
        parameters = self.request.body_arguments
        # Save filename into the input dictionary.
        parameters['input'] = filen
        
        if uptasks.allowed_file(filen) == 'True':
            with open(os.path.join(settingup.UPLOAD_FOLDER, filen), 'wb') as upload_file:
                upload_file.write(fileinfo['body'])
            
            pool = ThreadPoolExecutor(max_workers=settingmain.MAX_WORKERS)
            # if MySQL option is checked, execute the MySQL query script. Else, use the XML parsinh function.
            if 'mysql' in parameters:
                result = yield pool.submit(uniparser.uniprotParserSQL, parameters)
                pool.shutdown()
            else:
                result = yield pool.submit(uniparser.uniprot_parser, parameters)
                pool.shutdown()
            
        self.render("uniparser2_result.html", title="Uniprot Parser Result", result=result,)
        
settings = {
    "autoreload": True,
    "debug": True,
    "static_path": settingmain.APP_STATIC,
    "template_path": settingmain.APP_TEMPLATE,
    }

if __name__ == "__main__":
    application = tornado.web.Application([
        (r"/", MainHandler),
        (r"/ws", wsUI),
        (r"/websocket", WebSocketHandler),
        (r"/sequonparser", SequonParserUploadForm),
        (r"/sequonparser/process", SequonParserUploadFileProcess),
        (r"/sequonparser/results/download/(.*)", tornado.web.StaticFileHandler, {"path": "./sequon/results"}),
        (r"/uniprotparser", UniprotUploadForm),
        (r"/uniprotparser/process", UniprotUploadFileProcess),
        (r"/uniparser", UniparserForm),
        (r"/uniparser/process", UniparserProcess),
        (r"/uniprotparser/results/download/(.*)", tornado.web.StaticFileHandler, {"path": "./uniprot_parser/results"}),
        (r"/compare", CSVCompareForm),
        (r"/compare/process", CSVCompareProcess),
        (r"/compare/results/download/(.*)", tornado.web.StaticFileHandler, {"path": "./csvcompare/results"}),
        (r"/go", GOEnrichmentForm),
        (r"/go/process", GOEnrichmentResult),
        (r"/go/results/download/(.*)", tornado.web.StaticFileHandler, {"path": "./goenrichment/results"}),
        (r"/job_result/(?P<channel_id>[^\/]+)/(?P<job_key>[^\/]+)?", GetResult),
        (r"/sca", SubcellLocAnalyze),
        (r"/sca/process", SubcellLocAnalyzeProcess),
    ], **settings)
    application.listen(8888)
    tornado.ioloop.IOLoop.current().start()
