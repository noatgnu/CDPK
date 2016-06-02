from lxml import etree
#from xml.etree import ElementTree as etree
import os
import gzip
import csv
import re
from . import settingup
import MySQLdb
import time
from pyteomics import mass, parser
def mysql_cred():
    dbuser = settingup.DATABASES['default']['USER']
    dbpass = str(settingup.DATABASES['default']['PASSWORD'])
    dbhost = settingup.DATABASES['default']['HOST']
    db = settingup.DATABASES['default']['DB']
    return {'user':dbuser, 'pass':dbpass, 'host':dbhost, 'db':db}
def allowed_file(filename):
    # This function checked the extension of the input file against a list of allowed extensions.
    # Return True if it was in the list and False if it was not in the list.
    split_fn = filename.split('.')[-1] # Split filename at the '.' characters and only obtain the final item of the array.
    if split_fn:
       if split_fn in settingup.ALLOWED_EXTENSIONS:
           return 'True'
       else:
           return 'False'
    return 'False'

def check_column(input_path):
    # The program tried to discern the column which contain the uniprot accession id. 
    condition = 0
    accession_col_name = ''
    with open(input_path, 'rt') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel', delimiter='\t')
        for row in reader:
            for k in row:
                if re.search('\|[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\|', str(row[k])):
                    # First checking only the first row using regular expression for column containing the text with the id surrounded by 2 '|' separators.
                    accession_col_name = k
                    condition = 1
            if condition == 0:
                if re.search('^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', str(row[k])):
                    # If the first operation was not true, it checked for column with only the id and no separators.
                    accession_col_name = k
                    condition = 1
            break
    return accession_col_name

def uniprotParserSQL(parameters):
    output = dict()
    list_entry = set()
    no_match = []
    row = dict()
    colm = check_column(os.path.join(settingup.UPLOAD_FOLDER, parameters['input']))
    print(time.time())
    with open(os.path.join(settingup.UPLOAD_FOLDER, parameters['input']), 'rt') as infile:
        reader = csv.DictReader(infile, dialect='excel', delimiter='\t')
        for row in reader:
            result = re.search('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', str(row[colm]))
            list_entry.add(result.group(0))
    
    dbuser = settingup.DATABASES['default']['USER']
    dbpass = str(settingup.DATABASES['default']['PASSWORD'])
    dbhost = settingup.DATABASES['default']['HOST']
    db = settingup.DATABASES['default']['DB']
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=db)
    cursor = dbcon.cursor()
    db_choice = ''
    fieldname = ['accession']
    if 'dtype' in parameters:
        db_choice = settingup.DBTYPEMSQL[parameters['dtype'][0].decode('utf-8')]
    if 'database' in parameters:
        db = parameters['database'][0].decode('utf-8')
    if 'entryname' in parameters:
        fieldname.append('entry name')
    if 'proteinname' in parameters:
        fieldname.append('protein name')
    if 'subcellularlocation' in parameters:
        fieldname.append('subcellular location')
    if 'goid' in parameters:
        fieldname.append('geneontology IDs')
    if 'goterm' in parameters:
        fieldname.append('geneontology terms')
    if 'sequence' in parameters:
        fieldname.append('sequence')
    if 'sequencefeatures' in parameters:
        fieldname.append('sequence features')
    if 'glycosylation' in parameters:
        fieldname.append('glycosylation sites')
    
    gs_file = ''
    
    if 'gostats' in parameters:
        output['GOSTATS'] = 'gs_'+parameters['input']
        gs_file = open(os.path.join(settingup.RESULTS_FOLDER, 'gs_'+parameters['input']), 'wt')
        gs_empty = open(os.path.join(settingup.RESULTS_FOLDER, 'gs_empty_'+parameters['input']), 'wt')
        gswriter = csv.DictWriter(gs_file, fieldnames=['Geneontology IDs','Gocode','Entry'], dialect='excel', delimiter='\t')
        gswriter.writeheader()
        gs_e_writer = csv.DictWriter(gs_empty, fieldnames=['Geneontology IDs','Gocode','Entry'], dialect='excel', delimiter='\t')
        gs_e_writer.writeheader()
    pque = """SELECT * FROM """
    accq = """ WHERE accession = %s;"""
    restq = """ WHERE entryName = %s;"""
    with open(os.path.join(settingup.RESULTS_FOLDER, 'matched_'+parameters['input']), 'wt') as outfile:
        output['Main Output'] = 'matched_'+parameters['input']
        owriter = csv.DictWriter(outfile, fieldnames=fieldname, dialect='excel', delimiter='\t')
        owriter.writeheader()
        row = dict()
        for i in list_entry:
            
            for k in fieldname:
                row[k] = ''
            row['accession'] = i
            cursor.execute(pque+db+'Uniprot'+db_choice+'ACC'+accq, (i,))
            if not cursor.rowcount:
                no_match.append(i)
                continue
            else:
                result = cursor.fetchone()
                row['entry name'] = result[0]
            
            if 'proteinname' in parameters:
                cursor.execute(pque+db+'Uniprot'+db_choice+'FullName'+restq, (row['entry name'],))
                if cursor.rowcount:
                    result = cursor.fetchone()
                    row['protein name'] = result[1]
                cursor.execute(pque+db+'Uniprot'+db_choice+'EC'+restq, (row['entry name'],))
                if cursor.rowcount:
                    result = cursor.fetchone()
                    row['protein name'] += '('+result[1]+')'
            if 'subcellularlocation' in parameters:
                cursor.execute(pque+db+'Uniprot'+db_choice+'SubLoc'+restq, (row['entry name'],))
                if cursor.rowcount:
                    results = cursor.fetchall()
                    for r in results:
                        row['subcellular location'] += r[1]+':'+r[2]+';'
                        
            if 'goid' in parameters or 'goterm' in parameters:
                gostats = {'Entry':row['entry name'], 'Gocode':'IEA', 'Geneontology IDs':''}
                cursor.execute(pque+db+'Uniprot'+db_choice+'GO'+restq, (row['entry name'],))
                if cursor.rowcount:
                    results = cursor.fetchall()
                    for r in results:
                        if 'gostats' in parameters:
                            gostats['Geneontology IDs'] = r[1]
                            gswriter.writerow(gostats)
                        if 'goid' in parameters:
                            row['geneontology IDs'] += r[1]+';'
                        if 'goterm' in parameters:
                            row['geneontology terms'] += r[1]+'('+r[3]+');'
                else:
                    output['GOSTATS Empty'] = 'gs_empty_'+parameters['input']
                    if 'gostats' in parameters:
                        gs_e_writer.writerow(gostats)
            if 'sequence' in parameters:
                cursor.execute(pque+db+'Uniprot'+db_choice+'Seq'+restq, (row['entry name'],))
                if cursor.rowcount:
                    results = cursor.fetchone()
                    row['sequence'] = results[1]
            if 'sequencefeatures' in parameters or 'glycosylation' in parameters:
                cursor.execute(pque+db+'Uniprot'+db_choice+'FT'+restq, (row['entry name'],))
                if cursor.rowcount:
                    results = cursor.fetchall()
                    if 'sequencefeatures' in parameters:
                        for r in results:
                            row['sequence features'] += r[1]+':'+r[3]+'('+r[2]+');'
                    if 'glycosylation' in parameters:
                        for r in results:
                            if r[1] == 'glycosylation site':
                                row['glycosylation sites'] += r[1]+':'+r[3]+'('+r[2]+');'
            owriter.writerow(row)
    if 'gostats' in parameters:
        gs_file.close()
        gs_empty.close()
    if len(no_match) >= 0:
        with open(os.path.join(settingup.RESULTS_FOLDER, 'no_match_'+parameters['input']), 'wt') as no_match_file:
            output['No match'] = 'no_match_'+parameters['input']
            no_match_writer = csv.DictWriter(no_match_file, fieldnames=['accession'], dialect='excel', delimiter='\t')
            no_match_writer.writeheader()
            no_match_row = {'accession':''}
            for i in no_match:
                no_match_row['accession'] = i
                no_match_writer.writerow(no_match_row)
    print(time.time())
    return output
    
                        
    
    
def uniprot_parser(parameters):
    list_entry = set()
    
    colm = check_column(os.path.join(settingup.UPLOAD_FOLDER, parameters['input']))
    with open(os.path.join(settingup.UPLOAD_FOLDER, parameters['input']), 'rt') as infile:
        reader = csv.DictReader(infile, dialect='excel', delimiter='\t')
        for row in reader:
            result = re.search('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', str(row[colm]))
            list_entry.add(result.group(0))
    print(time.time())
    ns = '{http://uniprot.org/uniprot}'
    sub_list = ''
    database_path = ''
    db_choice = ''
    if 'dtype' in parameters:
        db_choice = settingup.DBTYPEMSQL[parameters['dtype'][0].decode('utf-8')]
    if 'database' in parameters:
        db = settingup.DBXML[parameters['database'][0].decode('utf-8')]
        database_path = os.path.join(settingup.DB_FOLDER, db_choice+'_'+db)
    go_found = False
    found_entry = False
    match_found = False
    relative_pos_acc = False
    relative_pos_name = False
    sub_loc = False
    recomname_loc = False
    goid_loc = False
    feature_loc = False
    feature_loc_seq = False
    dbr_id = False
    sub_listd = False
    go_idd = False
    go_termd = False
    seq_fd = False
    
    output = dict()
    protein = ''
    recommendedn = ''
    subcellular_loc = ''
    dbreference = ''
    feature_inf = ''
    position_seq = ''
    fieldname = ['accession']
    if 'entryname' in parameters:
        fieldname.append('entry name')
    if 'proteinname' in parameters:
        fieldname.append('protein name')
    if 'subcellularlocation' in parameters:
        fieldname.append('subcellular location')
    if 'goid' in parameters:
        fieldname.append('geneontology IDs')
    if 'goterm' in parameters:
        fieldname.append('geneontology terms')
    if 'sequence' in parameters:
        fieldname.append('sequence')
    if 'sequencefeatures' in parameters:
        fieldname.append('sequence features')
    if 'glycosylation' in parameters:
        fieldname.append('glycosylation sites')
    for k in fieldname:
        row[k] = ''
    with open(os.path.join(settingup.RESULTS_FOLDER, 'matched_'+parameters['input']), 'wt') as outfile:
        #fieldname = ['accession','entry name','protein name','subcellular location']
        output['Main Output'] = parameters['input']
        owriter = csv.DictWriter(outfile, fieldnames=fieldname, dialect='excel', delimiter='\t')
        owriter.writeheader()
        row = dict()   
        
        goid_list = ''
        goterm_list = ''
        seq_features = ''
        seq_position = ''
        ind_feature = ''
        gs_file = ''
        gswriter = ''
        gsrow = dict()
        goframe = ''
        gs_e = dict()
        
        if 'gostats' in parameters:
            gs_file = open(os.path.join(settingup.RESULTS_FOLDER, 'gs_'+parameters['input']), 'wt')
            gs_empty = open(os.path.join(settingup.RESULTS_FOLDER, 'gs_empty_'+parameters['input']), 'wt')
            
            
            gswriter = csv.DictWriter(gs_file, fieldnames=['Entry','Gocode','Geneontology IDs'], dialect='excel', delimiter='\t')
            gswriter.writeheader()
            gs_e_writer = csv.DictWriter(gs_empty, fieldnames=['Entry','Gocode','Geneontology IDs'], dialect='excel', delimiter='\t')
            gs_e_writer.writeheader()
                    
        for event, elem in etree.iterparse(gzip.open(database_path, 'rb'), events = ("start","end")):
            #print(event, elem.tag)
            
            if elem.tag == ns+'entry' and event == "start":
                found_entry = True 
                
            if elem.tag == ns+'accession' and elem.text in list_entry and found_entry == True and event == "end":
                match_found = True
                relative_pos_name = True
                row['accession'] = elem.text
                #print(elem.text)
                
                
            if elem.tag == ns+'name' and found_entry == True and match_found == True and relative_pos_name == True and event == "end" and 'entryname' in parameters:
                row['entry name'] = elem.text
                relative_pos_name = False
                
            if elem.tag == ns+'recommendedName' and match_found == True and found_entry == True and event == "start" and 'proteinname' in parameters:
                recomname_loc = True
                
            if elem.tag == ns+'fullName' and match_found == True and recomname_loc == True and event == "end" and 'proteinname' in parameters:
                row['protein name'] = elem.text
                
            if elem.tag == ns+'ecNumber' and match_found == True and recomname_loc == True and event == "end" and 'proteinname' in parameters:
                row['protein name'] += ' (EC '+ elem.text + ')'
                
            if elem.tag == ns+'recommendedName' and match_found == True and found_entry == True and event == "end" and 'proteinname' in parameters:
                recomname_loc = False
                
            if elem.tag == ns+'comment' and match_found == True and elem.get('type') == 'subcellular location' and found_entry == True and event == "start" and 'subcellularlocation' in parameters:
                sub_loc = True
                #print(etree.tostring(elem,pretty_print=True))
                
            if elem.tag == ns+'location' and match_found == True and found_entry == True and sub_loc == True and event == "end" and 'subcellularlocation' in parameters:
                if sub_listd == True:
                    sub_list = sub_list+';'
                sub_list = sub_list+elem.text
                sub_listd = True
            if elem.tag == ns+'topology' and match_found == True and found_entry == True and sub_loc == True and event == "end" and 'subcellularlocation' in parameters:
                if sub_listd == True:
                    sub_list = sub_list+';'
                sub_list = sub_list+elem.text
                sub_listd = True
            if elem.tag == ns+'text' and match_found == True and found_entry == True and sub_loc == True and event == "end" and 'subcellularlocation' in parameters:
                if sub_listd == True:
                    sub_list = sub_list+';'
                sub_list = sub_list+elem.text
                sub_listd = True
            
            if elem.tag == ns+'comment' and match_found == True and elem.get('type') == 'subcellular location' and found_entry == True and event == "end" and 'subcellularlocation' in parameters:
                sub_loc = False
                row['subcellular location'] = sub_list
                sub_list = ''
                sub_listd = False
                #print(etree.tostring(elem, pretty_print=True))
            if elem.tag == ns+'dbReference' and match_found == True and found_entry == True and elem.get('type') == 'GO' and event == "start":
                go_found = True
                #print(go_idd)
                if 'goid' in parameters or 'gostats' in parameters:
                    goid_loc = True
                    if go_idd == True:
                        goid_list = goid_list + ';' + elem.get('id')
                    else:
                        goid_list = elem.get('id')
                        go_idd = True
                     
                if 'goterm' in parameters:
                    goframe = elem.get('id')
                
            if elem.tag == ns+'property' and match_found == True and found_entry == True and goid_loc == True and 'goterm' in parameters and event == "end":
                if elem.get('type') == 'term':
                    if go_termd == True:
                        goterm_list = goterm_list+';'
                        
                    goterm_list = goterm_list+goframe+'('+elem.get('value')+')'
                    go_termd = True
            if elem.tag == ns+'dbReference' and match_found == True and found_entry == True and elem.get('type') == 'GO' and event == "end":
                if 'goid' in parameters or 'gostats' in parameters:
                    goid_loc = False
                
            if elem.tag == ns+'feature' and match_found == True and found_entry == True and event == "start":
                if 'sequencefeatures' in parameters or 'glycosylation' in parameters:
                    feature_loc = True
                    if seq_fd == True:
                        seq_features = seq_features + ';'
                    seq_fd = True
                    if 'description' in elem.attrib:
                        ind_feature = elem.get('type')+ ': ' + elem.get('description')
                    else:
                        ind_feature = elem.get('type')
                    
                    
            if elem.tag == ns+'location' and match_found == True and found_entry == True and event == "start":
                if 'goid' in parameters or 'gostats' in parameters:
                    feature_loc_seq = True
                
            if elem.tag == ns+'begin' and match_found == True and found_entry == True and feature_loc_seq == True and event == "end" and 'sequencefeatures' in parameters:
                
                if not 'position' in elem.attrib:
                    seq_position = 'unknown'
                        
                else:
                    seq_position = elem.get('position')
                                
            if elem.tag == ns+'end' and match_found == True and found_entry == True and feature_loc_seq == True and event == "end" and 'sequencefeatures' in parameters:
                
                if not 'position' in elem.attrib:
                    seq_position = seq_position +'-'+ 'unknown'
                        
                else:
                    seq_position = seq_position +'-'+elem.get('position')
                
            if elem.tag == ns+'position' and match_found == True and found_entry == True and feature_loc_seq == True and event == "end":
                if 'sequencefeatures' in parameters or 'glycosylation' in parameters:
                    if not 'position' in elem.attrib:
                        seq_position = 'unknown'
                            
                    else:
                        seq_position = elem.get('position')
                    
            if elem.tag == ns+'location' and match_found == True and found_entry == True and event == "end":
                if 'sequencefeatures' in parameters or 'glycosylation' in parameters:
                    feature_loc_seq = False
                
            if elem.tag == ns+'feature' and match_found == True and found_entry == True and event == "end":
                if 'goid' in parameters or 'gostats' in parameters:
                    feature_loc = False
                    seq_features = seq_features + ind_feature+'('+seq_position+')'
                
            if elem.tag == ns+'sequence' and match_found == True and found_entry == True and event == "end" and 'sequence' in parameters:
                row['sequence'] = elem.text
                
            if elem.tag == ns+'entry' and match_found == True and found_entry == True and event == "end":
                
                #print(goid_list)
                if go_found == True:
                    output['GOSTATS'] = 'gs_'+parameters['input']
                    if 'goid' in parameters:
                        row['geneontology IDs'] = goid_list
                    if 'goterm' in parameters:
                        row['geneontology terms'] = goterm_list
                else:
                    output['GOSTATS Entry without any GOIDs'] = 'gs_empty_'+parameters['input']
                    if 'goid' in parameters:
                        row['geneontology IDs'] = ''
                    if 'goterm' in parameters:
                        row['geneontology terms'] = ''
                    gs_e['Entry'] = row['accession']
                    gs_e['Gocode'] = 'IEA'
                    gs_e['Geneontology IDs'] = ''
                    gs_e_writer.writerow(gs_e)
                    
                    
                if 'sequencefeatures' in parameters:
                    row['sequence features'] = seq_features
                if 'glycosylation' in parameters:
                    row['glycosylation sites'] = seq_features
                if 'gostats' in parameters:
                    if  'entry name' in row:
                        if not row['entry name'] == '':
                            gsrow['Entry'] = row['entry name']
                        else:
                            gsrow['Entry'] = row['accession']
                    else:
                        gsrow['Entry'] = row['accession']
                    goid_ar = goid_list.split(';')
                    #print(goid_ar)
                    for goid in goid_ar:
                        gsrow['Geneontology IDs'] = goid
                        gsrow['Gocode'] = 'IEA'
                        #print(gsrow)
                        if go_found == True:
                            gswriter.writerow(gsrow)
                owriter.writerow(row)
                list_entry.remove(row['accession'])    
                #print(len(list_entry))
               
                goid_list = ''
                goterm_list = ''
                seq_features = ''
                go_idd = False
                go_termd = False
                seq_fd = False
                gs_e = dict()
                gsrow = dict()
                elem.clear()
                go_found = False
                found_entry = False
                match_found = False
                for k in fieldname:
                    row[k] = ''
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
                if len(list_entry) == 0:
                    if 'gostats' in parameters:
                        gs_file.close()
                        gs_empty.close()
                    print(time.time())
                    return output                    
            if elem.tag == ns+'entry' and found_entry == False and event == "end":
                found_entry = False
                
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
        if 'gostats' in parameters:
            gs_file.close()
            gs_empty.close()
        if not len(list_entry) == 0:
            output['No Match'] = 'no_match_'+parameters['input']
            field_mini = ['accession']
            with open(os.path.join(settingup.RESULTS_FOLDER, 'no_match_'+parameters['input']), 'wt') as gs_file_no_hit:
                no_match = csv.DictWriter(gs_file_no_hit, fieldnames=field_mini, dialect='excel', delimiter='\t')
                no_match.writeheader()
                no_match_row = dict()
                for i in list_entry:
                    no_match_row['accession'] = i
                    no_match.writerow(no_match_row)
    print(time.time())
    return output
    
def compdb(filepath='uniprot_sprot.xml.gz', output=settingup.DBXML['Compact']):    
    with gzip.open(output, 'wb', 9) as outfile:
        found_entry = False
        relative_pos_acc = False
        relative_pos_name = False
        sub_loc = False
        recomname_loc = False
        goid_loc = False
        feature_loc = False
        feature_loc_seq = False
        dbr_id = False
        ns = '{http://uniprot.org/uniprot}'
        entry = etree.Element('entry')
        protein = ''
        recommendedn = ''
        subcellular_loc = ''
        #subr_loc = etree.Element('subcellularlocation')
        dbreference = ''
        feature_inf = ''
        position_seq = ''
        found_organism = False
        #fieldname = ['accession','entry name','protein name','subcellular location']
        #owriter = csv.DictWriter(outfile, fieldnames=fieldname, dialect='excel', delimiter='\t')
        #owriter.writeheader()
        outfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        outfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
        #row = dict()    
        for event, elem in etree.iterparse(gzip.open(filepath, 'rb'), events = ("start","end")):
            #print(event, etree.tostring(elem))
            #if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and event == 'start':
                #print(etree.tostring(elem))
            
            if elem.tag == ns+'entry' and event == "start":
                found_entry = True 
                
            if elem.tag == ns+'accession' and found_entry == True and event == "end":
                
                accession = etree.SubElement(entry, 'accession')
                accession.text = elem.text
                #print(elem.text)
                #print(etree.tostring(entry))
                relative_pos_name = True
            if elem.tag == ns+'name' and found_entry == True and relative_pos_name == True and event == "end":
                e_name = etree.SubElement(entry,'name')
                e_name.text = elem.text
                
                relative_pos_name = False
            if elem.tag == ns+'protein' and found_entry == True and event == "start":
                protein = etree.SubElement(entry,'protein')
            if elem.tag == ns+'recommendedName' and found_entry == True and event == "start":
                recomname_loc = True
                recommendedn = etree.SubElement(protein, 'recommendedName')
            if elem.tag == ns+'fullName' and recomname_loc == True and event == "end":
                fullname = etree.SubElement(recommendedn, 'fullName')
                fullname.text = elem.text
            if elem.tag == ns+'ecNumber' and recomname_loc == True and event == "end":
                ecnm = etree.SubElement(recommendedn, 'ecNumber')
                ecnm.text = elem.text
            if elem.tag == ns+'recommendedName' and found_entry == True and event == "end":
                recomname_loc = False
            if elem.tag == ns+'organism' and found_entry == True and event == "start":
                found_organism = True
                organism = etree.SubElement(entry, 'organism')
            if elem.tag == ns+'name' and found_entry == True and elem.get('type') == 'scientific' and found_organism == True and event == "end":
                org_name = etree.SubElement(organism, 'name', {'type':'scientific'})
                org_name.text = elem.text
            if elem.tag == ns+'organism' and found_entry == True and event == "end":
                found_organism = False
            if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "start":
                sub_loc = True
                subcellular_loc = etree.SubElement(entry, 'comment', {'type': 'subcellular location'})
            if elem.tag == ns+'subcellularLocation' and found_entry == True and event == "start":
                
                subr_loc = etree.SubElement(subcellular_loc, 'subcellularLocation')
                  
            if elem.tag == ns+'location' and found_entry == True and sub_loc == True and event == "end":
                location = etree.SubElement(subr_loc, 'location')
                location.text = elem.text
                #print(elem.text)
            if elem.tag == ns+'topology' and found_entry == True and sub_loc == True and event == "end":
                topology = etree.SubElement(subr_loc, 'topology')
                topology.text = elem.text
            if elem.tag == ns+'orientation' and found_entry == True and sub_loc == True and event == "end":
                
                topology = etree.SubElement(subr_loc, 'orientation')
                topology.text = elem.text
                
            if elem.tag == ns+'text' and found_entry == True and sub_loc == True and event == "end":
                text = etree.SubElement(subcellular_loc, 'text')
                text.text = elem.text
                
                #print(elem.text)
            #if elem.tag == ns+'subcellular location' and found_entry == True and event == "end":
                #subcellular_loc.append(subr_loc)
                #subr_loc = etree.Element('subcellular location')
                
            if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "end":
                sub_loc = False
                
            if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "start":
                goid_loc = True
                dbre_at = {'type':elem.get('type'), 'id':elem.get('id')}
                dbReference = etree.SubElement(entry, 'dbReference', dbre_at)
                
            if elem.tag == ns+'property' and found_entry == True and goid_loc == True and event == "end":
                dbre_at = {'type':elem.get('type'), 'value':elem.get('value')}
                dbrprop = etree.SubElement(dbReference, 'property', dbre_at)
            if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "end":
                goid_loc = False
                
            if elem.tag == ns+'feature' and found_entry == True and event == "start":
                feature_loc = True
                fattribute = {'type':elem.get('type')}
                if elem.get('description'):
                    fattribute['description'] = elem.get('description')
                    
                feature_inf = etree.SubElement(entry, 'feature', fattribute)
                
                
            if elem.tag == ns+'location' and found_entry == True and event == "start":
                position_seq = etree.SubElement(feature_inf, 'location')
                feature_loc_seq = True
            if elem.tag == ns+'begin' and found_entry == True and feature_loc_seq == True and event == "end":
                b_start = dict()
                if not 'position' in elem.attrib:
                    b_start['status'] = 'unknown'
                        
                else:
                    b_start['position'] = elem.get('position')
                
                b_pos = etree.SubElement(position_seq, 'begin', b_start)
                
            if elem.tag == ns+'end' and found_entry == True and feature_loc_seq == True and event == "end":
                b_end = dict()
                if not 'position' in elem.attrib:
                    b_end['status'] = 'unknown'
                        
                else:
                    b_end['position'] = elem.get('position')
                e_pos = etree.SubElement(position_seq, 'end', b_end)
                
            if elem.tag == ns+'position' and found_entry == True and feature_loc_seq == True and event == "end":
                b_direct = dict()
                if not 'position' in elem.attrib:
                    b_direct['status'] = 'unknown'
                        
                else:
                    b_direct['position'] = elem.get('position')
                pos = etree.SubElement(position_seq, 'position', b_direct)
                
            if elem.tag == ns+'location' and found_entry == True and event == "end":
                feature_loc_seq = False
                
            if elem.tag == ns+'feature' and found_entry == True and event == "end":
                feature_loc = False
            if elem.tag == ns+'sequence' and found_entry == True and event == "end" and 'checksum' in elem.attrib:
                sattribute = dict()
                if elem.get('length'):
                    sattribute['length'] = elem.get('length')
                if elem.get('mass'):
                    sattribute['mass'] = elem.get('mass')
                sequence = etree.SubElement(entry, 'sequence', sattribute)
                
                sequence.text = elem.text
            if elem.tag == ns+'entry' and found_entry == True and event == "end":
                found_entry = False
                
                #print(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                outfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                elem.clear()
                entry = etree.Element('entry')
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
            if elem.tag == ns+'entry' and found_entry == False and event == "end":
                found_entry = False
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
        outfile.write(b'</uniprot>')
 
def compdb_sql(filepath='uniprot_sprot.xml.gz', dtype='sp', compact=True, compact_only=False, taxa_list=settingup.TAXA):
    count = 0
    p1n = 'Uniprot'
    taxa_dict = dict()
    cn = 'Compact'
    ac1 = p1n+dtype+'ACC'
    fn1 = p1n+dtype+'FullName'
    ec1 = p1n+dtype+'EC'
    or1 = p1n+dtype+'Organism'
    sl1 = p1n+dtype+'SubLoc'
    go1 = p1n+dtype+'GO'
    ft1 = p1n+dtype+'FT'
    se1 = p1n+dtype+'Seq'
    tables = []
    tables = taxa_list
    for i in taxa_list:
        taxa_dict[i] = False
    
    if compact_only == True:
        tables = [cn]
    if compact == True:
        tables.append(cn)
        
    dbuser = settingup.DATABASES['default']['USER']
    dbpass = str(settingup.DATABASES['default']['PASSWORD'])
    dbhost = settingup.DATABASES['default']['HOST']
    db = settingup.DATABASES['default']['DB']
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=db)
    cursor = dbcon.cursor()
    for tax in tables:
        dacc = """DROP TABLE IF EXISTS """ + tax+ac1
        cursor.execute(dacc)
        dfn = """DROP TABLE IF EXISTS """ + tax+fn1
        cursor.execute(dfn)
        dec = """DROP TABLE IF EXISTS """ + tax+ec1
        cursor.execute(dec)
        dorg = """DROP TABLE IF EXISTS """ + tax+or1
        cursor.execute(dorg)
        dsubl = """DROP TABLE IF EXISTS """ + tax+sl1
        cursor.execute(dsubl)
        dgo = """DROP TABLE IF EXISTS """ + tax+go1
        cursor.execute(dgo)
        dft = """DROP TABLE IF EXISTS """ + tax+ft1
        cursor.execute(dft)
        dseq = """DROP TABLE IF EXISTS """ + tax+se1
        cursor.execute(dseq)
        tbc_acc = """CREATE TABLE IF NOT EXISTS """+ tax+ac1 +""" (entryName char(20) NOT NULL, accession char(20) NOT NULL);""" 
        cursor.execute(tbc_acc)
        tbc_fn = """CREATE TABLE IF NOT EXISTS """+ tax+fn1 +""" (entryName char(20) NOT NULL, fullName text NOT NULL);""" 
        cursor.execute(tbc_fn)
        tbc_ec = """CREATE TABLE IF NOT EXISTS """+tax+ec1+""" (entryName char(20) NOT NULL, EC char(20) NOT NULL);""" 
        cursor.execute(tbc_ec)
        tbc_org = """CREATE TABLE IF NOT EXISTS """+tax+or1+""" (entryName char(20) NOT NULL, organismName text NOT NULL);""" 
        cursor.execute(tbc_org)
        tbc_sub_loc = """CREATE TABLE IF NOT EXISTS """+tax+sl1+""" (entryName char(20) NOT NULL, subType char(20) NOT NULL, location text NOT NULL);"""
        cursor.execute(tbc_sub_loc)
        tbc_go = """CREATE TABLE IF NOT EXISTS """+tax+go1+""" (entryName char(20) NOT NULL, GOID char(20) NOT NULL, evidence char(255) NOT NULL, GOTerm text NOT NULL);"""
        cursor.execute(tbc_go)
        tbc_ft = """CREATE TABLE IF NOT EXISTS """+tax+ft1+""" (entryName char(20) NOT NULL, FTType text NOT NULL, position char(50) NOT NULL, description text NOT NULL);"""
        cursor.execute(tbc_ft)
        tbc_seq = """CREATE TABLE IF NOT EXISTS """+tax+se1+""" (entryName char(20) NOT NULL, sequence text NOT NULL);"""
        cursor.execute(tbc_seq)
    
    ins_p = """INSERT INTO """
    up_p = """UPDATE """
    ins_acc = """ (entryName, accession) VALUES (%s, %s);"""
    ins_fn = """ (entryName, fullName) VALUES (%s, %s);"""
    ins_ec = """ (entryName, EC) VALUES (%s, %s);"""
    ins_org = """ (entryName, organismName) VALUES (%s, %s);"""
    ins_sub_loc = """ (entryName, subType, location) VALUES (%s, %s, %s);"""
    ins_go = """ (entryName, GOID, evidence, GOTerm) VALUES (%s, %s, %s, %s);"""
    ins_ft = """ (entryName, FTType, position, description) VALUES (%s, %s, %s, %s);"""
    ins_seq = """ (entryName, sequence) VALUES (%s, %s);"""
    
    qaccst = []
    qfnst = []
    qecst = []
    qorgst = []
    qsublst = []
    qgost = []
    qftst = []
    qseqst = []
    sname_found = False
    submitted_name = ''
    found_entry = False # Checkpoint for each uniprot entry.
    relative_pos_name = False # Checkpoint for the first recommended name's full name entry.
    sub_loc = False # Checkpoint for subcellular location.
    recomname_loc = False # Checkpoint for recommended name.
    goid_loc = False # Checkpoint for GOID and GOTerm.
    feature_loc = False # Checkpoint for feature entries.
    feature_loc_seq = False # Checkpoint for feature location.
    dbr_id = False # Checkpoint for DBreference entries.
    found_organism = False # Checkpoint for organism entries.
    ns = '{http://uniprot.org/uniprot}' # Setting up namespace for the Uniprot XML document.
    # Setting up a dictionary which would temporary hold for all the parsed information before write out.
    entry = {'entryName':'', 'accessions':[], 'ec':'', 'fullName':'', 'organismName':'', 'subcellLoc':[], 'go':[], 'ft':[], 'seq':''}

    # Setting up temporary holder for each subcellular location information entry.
    sl_stuff = [] #{'type':'', 'location':''}
    # Setting up temporary holder for each GOID and GOTerm information entry.
    go_stuff = [] #{'id':'', 'evidence':'', 'GOTerm':''}
    # Setting up temporary holder for each sequence feature.
    ft_stuff = [] #{'type':'', 'position':'', 'description':''}
    ftpos = ''
    # Begin parsing using lxml's iterparse.
    for event, elem in etree.iterparse(gzip.open(filepath, 'rb'), events = ("start","end")):
        
        if elem.tag == ns+'entry' and event == "start":
            # Start entry checkpoint.
            found_entry = True
            accession_list = []
            
        if elem.tag == ns+'accession' and found_entry == True and event == "end":
            entry['accessions'].append(elem.text)
            relative_pos_name = True
            # Start entry name checkpoint.
            # Add each found accession into the temporary entry dictionary.
            
        if elem.tag == ns+'name' and found_entry == True and relative_pos_name == True and event == "end":
            entry['entryName'] = elem.text
            relative_pos_name = False
            # End entry name checkpoint.
            # Add entry name to the temporary entry dictionary.
            
        if elem.tag == ns+'recommendedName' and found_entry == True and event == "start":
            recomname_loc = True
            # Start recommended name checkpoint.
            
        if elem.tag == ns+'fullName' and recomname_loc == True and event == "end":
            entry['fullName'] = elem.text
            recomname_loc = False
            
            # End recommeded name checkpoint.
            # Add full name into the temporary entry dictionary.
        if elem.tag == ns+'submittedName' and recomname_loc == False and found_entry == True and event == "start":
            sname_found = True
        if elem.tag == ns+'fullName' and sname_found == True and event == "end":
            submitted_name = elem.text
            sname_found = False
        if elem.tag == ns+'ecNumber' and found_entry == True and event == "end":
            entry['ec'] = elem.text
            # Add EC Number to the entry to the entry dictionary.
            
        if elem.tag == ns+'recommendedName' and found_entry == True and event == "end":
            recomname_loc = False
            # End recommeded name checkpoint.
            
        if elem.tag == ns+'organism' and found_entry == True and event == "start":
            found_organism = True
            # Start organism checkpoint.
            
        if elem.tag == ns+'name' and compact_only == False and found_entry == True and elem.get('type') == 'scientific' and found_organism == True and event == "end":
            entry['organismName'] = elem.text
            # Add organism name to the entry dictionary.
            
        if elem.tag == ns+'taxon' and found_organism == True and found_entry == True and event == "start" and elem.text in taxa_dict:
            taxa_dict[elem.text] = True
            # Turning on trigger for input into the appropriate taxon.
            
        if elem.tag == ns+'organism' and found_entry == True and event == "end":
            found_organism = False
            # End organism checkpoint.
            
        if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "start":
            sub_loc = True
            # Start subcellular location checkpoint.
            
        if elem.tag == ns+'location' and found_entry == True and sub_loc == True and event == "end":
            sl_stuff.append('location') # Assign subcellular location type to temporary subcell location dictionary.
            sl_stuff.append(elem.text) # Assign subcellular location description to temporary subcell location dictionary.
            #print(sl_stuff)
            entry['subcellLoc'].append(sl_stuff) # Send temporary subcell location dictionary to the entry temporary dictionary. 
            #print(entry)
            sl_stuff = [] # Wipe temporary subcell location dictionary.
            
        if elem.tag == ns+'topology' and found_entry == True and sub_loc == True and event == "end":
            sl_stuff.append('topology')
            sl_stuff.append(elem.text)
            
            entry['subcellLoc'].append(sl_stuff)
            sl_stuff = []
        if elem.tag == ns+'orientation' and found_entry == True and sub_loc == True and event == "end":
            sl_stuff.append('orientation')
            sl_stuff.append(elem.text)
            entry['subcellLoc'].append(sl_stuff)
            sl_stuff = []
        if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "end":
            sub_loc = False
            
        if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "start":
            
            goid_loc = True
            go_stuff = []
            if 'id' in elem.attrib:
                go_stuff.append(elem.get('id'))
                #print(elem.get('id'))
            
        if elem.tag == ns+'property' and found_entry == True and goid_loc == True and event == "end":
            if elem.get('type') == 'term':
                go_stuff.append(elem.get('value'))
            
            if elem.get('type') == 'evidence':
                go_stuff.append(elem.get('value'))
            #print(go_stuff)
                
        if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "end":
            goid_loc = False
            #print(go_stuff)
            entry['go'].append(go_stuff)
                
        if elem.tag == ns+'feature' and found_entry == True and event == "start":
            feature_loc = True
            ft_stuff = []
            ft_stuff.append(elem.get('type'))
            if 'description' in elem.attrib:
                ft_stuff.append(elem.get('description'))
            else:
                ft_stuff.append('')
            
        if elem.tag == ns+'location' and found_entry == True and event == "start":
            feature_loc_seq = True
            
        if elem.tag == ns+'begin' and found_entry == True and feature_loc_seq == True and event == "end":
            if not 'position' in elem.attrib:
                ft_stuff.append('unknown')
                    
            else:
                ft_stuff.append(elem.get('position'))
            
        if elem.tag == ns+'end' and found_entry == True and feature_loc_seq == True and event == "end":
            if not 'position' in elem.attrib:
                ft_stuff[2] += '-unknown'
            else:
                ft_stuff[2] += '-'+ elem.get('position')
            
        if elem.tag == ns+'position' and found_entry == True and feature_loc_seq == True and event == "end":
            if not 'position' in elem.attrib:
                ft_stuff.append('unknown')
            else:
                ft_stuff.append(elem.get('position'))
            
        if elem.tag == ns+'location' and found_entry == True and event == "end":
            feature_loc_seq = False
            
        if elem.tag == ns+'feature' and found_entry == True and event == "end":
            feature_loc = False
            entry['ft'].append(ft_stuff)
            
        if elem.tag == ns+'sequence' and found_entry == True and event == "end" and 'checksum' in elem.attrib:
            entry['seq'] = elem.text
        if elem.tag == ns+'entry' and found_entry == True and event == "end":
            count+=1
            found_entry = False
            if entry['fullName'] == '':
                entry['fullName'] = submitted_name
            
            #print(entry)
            print(count, end='\r')
            
            if compact_only == True or compact == True:
                
                for acc in entry['accessions']:
                    #cursor.execute(ins_p+cn+ac1+ins_acc, (entry['entryName'], acc))
                    #dbcon.commit()
                    qaccst.append((entry['entryName'], acc))
                if len(qaccst) >= 5000:
                    cursor.executemany(ins_p+cn+ac1+ins_acc, qaccst)
                    dbcon.commit()
                    qaccst = []
                    
                #cursor.execute(ins_p+cn+fn1+ins_fn, (entry['entryName'], entry['fullName']))
                #dbcon.commit()
                qfnst.append((entry['entryName'], entry['fullName']))
                if len(qfnst) >= 5000:
                    cursor.executemany(ins_p+cn+fn1+ins_fn, qfnst)
                    dbcon.commit()
                    qfnst = []
                #cursor.execute(ins_p+cn+ec1+ins_ec, (entry['entryName'], entry['ec']))
                #dbcon.commit()
                qecst.append((entry['entryName'], entry['ec']))
                if len(qecst) >= 5000:
                    cursor.executemany(ins_p+cn+ec1+ins_ec, qecst)
                    dbcon.commit()
                    qecst = []
                #cursor.execute(ins_p+cn+or1+ins_org, (entry['entryName'], entry['organismName']))
                #dbcon.commit()
                qorgst.append((entry['entryName'], entry['organismName']))
                if len(qorgst) >= 5000:
                    cursor.executemany(ins_p+cn+or1+ins_org, qorgst)
                    dbcon.commit()
                    qorgst = []
                if not len(entry['subcellLoc']) == 0:
                    for s in entry['subcellLoc']:
                        #cursor.execute(ins_p+cn+sl1+ins_sub_loc, (entry['entryName'], s[0], s[1]))
                        #dbcon.commit()
                        qsublst.append((entry['entryName'], s[0], s[1]))
                if len(qsublst) >= 5000:        
                    cursor.executemany(ins_p+cn+sl1+ins_sub_loc, qsublst)
                    dbcon.commit()
                    qsublst = []
                    
                if not len(entry['go']) == 0:
                    for g in entry['go']:
                        #cursor.execute(ins_p+cn+go1+ins_go, (entry['entryName'], g[0], g[2], g[1]))
                        #dbcon.commit()
                        qgost.append((entry['entryName'], g[0], g[2], g[1]))
                if len(qgost) >= 5000:
                    cursor.executemany(ins_p+cn+go1+ins_go, qgost)
                    dbcon.commit()
                    qgost = []
                    
                if not len(entry['ft']) == 0:
                    for f in entry['ft']:
                        #cursor.execute(ins_p+cn+ft1+ins_ft, (entry['entryName'], f[0], f[2], f[1]))
                        #dbcon.commit()
                        qftst.append((entry['entryName'], f[0], f[2], f[1]))
                if len(qftst) >= 5000:
                    cursor.executemany(ins_p+cn+ft1+ins_ft, qftst)
                    dbcon.commit()
                    qftst = []
                #cursor.execute(ins_p+cn+se1+ins_seq, (entry['entryName'], entry['seq']))
                #dbcon.commit()
                qseqst.append((entry['entryName'], entry['seq']))
                if len(qseqst) >= 5000:
                    cursor.executemany(ins_p+cn+se1+ins_seq, qseqst)
                    dbcon.commit()
                    qseqst = []
                    
            if compact_only == False:
                for t in taxa_dict:
                    if taxa_dict[t] == True:
                        for acc in entry['accessions']:
                            cursor.execute(ins_p+t+ac1+ins_acc, (entry['entryName'], acc))
                            dbcon.commit()
                        cursor.execute(ins_p+t+fn1+ins_fn, (entry['entryName'], entry['fullName']))
                        dbcon.commit()
                        cursor.execute(ins_p+t+ec1+ins_ec, (entry['entryName'], entry['ec']))
                        dbcon.commit()
                        cursor.execute(ins_p+t+or1+ins_org, (entry['entryName'], entry['organismName']))
                        dbcon.commit()
                        if not len(entry['subcellLoc']) == 0:
                            for s in entry['subcellLoc']:
                                cursor.execute(ins_p+t+sl1+ins_sub_loc, (entry['entryName'], s[0], s[1]))
                                dbcon.commit()
                        if not len(entry['go']) == 0:
                            for g in entry['go']:
                                cursor.execute(ins_p+t+go1+ins_go, (entry['entryName'], g[0], g[2], g[1]))
                                dbcon.commit()
                        if not len(entry['ft']) == 0:
                            for f in entry['ft']:
                                cursor.execute(ins_p+t+ft1+ins_ft, (entry['entryName'], f[0], f[2], f[1]))
                                dbcon.commit()
                        cursor.execute(ins_p+t+se1+ins_seq, (entry['entryName'], entry['seq']))
                        dbcon.commit()
            elem.clear()
            entry = {'entryName':'', 'accessions':[], 'ec':'', 'fullName':'', 'organismName':'', 'subcellLoc':[], 'go':[], 'ft':[], 'seq':''}
            submitted_name = ''
            while elem.getprevious() is not None:
                del elem.getparent()[0]
        if elem.tag == ns+'entry' and found_entry == False and event == "end":
            found_entry = False
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
    if not qaccst == []:
        cursor.executemany(ins_p+cn+ac1+ins_acc, qaccst)
        dbcon.commit()
        qaccst = []
    if not qfnst == []:
        cursor.executemany(ins_p+cn+fn1+ins_fn, qfnst)
        dbcon.commit()
        qfnst = []
    if not qecst == []:
        cursor.executemany(ins_p+cn+ec1+ins_ec, qecst)
        dbcon.commit()
        qecst = []
    if not qorgst == []:
        cursor.executemany(ins_p+cn+or1+ins_org, qorgst)
        dbcon.commit()
        qorgst = []
    if not qsublst == []:
        cursor.executemany(ins_p+cn+sl1+ins_sub_loc, qsublst)
        dbcon.commit()
        qsublst = []
    if not qgost == []:
        cursor.executemany(ins_p+cn+go1+ins_go, qgost)
        dbcon.commit()
        qgost = []
    if not qftst == []:
        cursor.executemany(ins_p+cn+ft1+ins_ft, qftst)
        dbcon.commit()
        qftst = []
    if not qseqst == []:
        cursor.executemany(ins_p+cn+se1+ins_seq, qseqst)
        dbcon.commit()
        qseqst = []
    if compact_only == False:
        mysqlindex_table(tables, dtype)
    if compact_only == True:
        mysqlindex_table(['Compact'], dtype)
    cursor.close()
    dbcon.close()
    print(count)
def mysqlindex_table(tables, dtype):
    dbuser = settingup.DATABASES['default']['USER']
    dbpass = str(settingup.DATABASES['default']['PASSWORD'])
    dbhost = settingup.DATABASES['default']['HOST']
    db = settingup.DATABASES['default']['DB']
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=db)
    cursor = dbcon.cursor()
    
    for t in tables:
        cursor.execute("""CREATE INDEX """+t+dtype+"""acc ON """+t+"""Uniprot"""+dtype+"""ACC (accession)""")
        dbcon.commit()
        cursor.execute("""CREATE INDEX """+t+dtype+"""fn ON """+t+"""Uniprot"""+dtype+"""FullName (entryName)""")
        dbcon.commit()
        cursor.execute("""CREATE INDEX """+t+dtype+"""ec ON """+t+"""Uniprot"""+dtype+"""EC (entryName)""")
        dbcon.commit()
        cursor.execute("""CREATE INDEX """+t+dtype+"""org ON """+t+"""Uniprot"""+dtype+"""Organism (entryName)""")
        dbcon.commit()
        cursor.execute("""CREATE INDEX """+t+dtype+"""sub ON """+t+"""Uniprot"""+dtype+"""SubLoc (entryName)""")
        dbcon.commit()
        cursor.execute("""CREATE INDEX """+t+dtype+"""go ON """+t+"""Uniprot"""+dtype+"""GO (entryName)""")
        dbcon.commit()
        cursor.execute("""CREATE INDEX """+t+dtype+"""ft ON """+t+"""Uniprot"""+dtype+"""FT (entryName)""")
        dbcon.commit()
        cursor.execute("""CREATE INDEX """+t+dtype+"""seq ON """+t+"""Uniprot"""+dtype+"""Seq (entryName)""")
        dbcon.commit()
    cursor.close()
    dbcon.close()
#def mysqlcleanup_list()
def compdb_taxa(filepath='uniprot_sprot.xml.gz', output=settingup.DBXML['Compact'], dtype='sp', taxa_list=['Vertebrata','Mammalia','Eukaryota','Fungi']):
    # Compacting and separating Uniprot XML database into separated taxa defined in "settingup.py" while also output a compact version of the full database. 
    # The compact version involved entry name, protein name, accessions, feature, sequence feature, goid and goterm, sequence.
    taxa_dict = dict()  
    for i in taxa_list: # taxa dictionary derives keys from the taxa_list array. Set default values of each key to false.
        taxa_dict[i] = False
    count = 0 # Start counting of entries.
    with gzip.open(os.path.join(settingup.DB_FOLDER, dtype+'_'+output), 'wb', 9) as outfile, gzip.open(os.path.join(settingup.DB_FOLDER, dtype+'_'+settingup.DBXML['Vertebrata']), 'wb', 9) as vboutfile, gzip.open(os.path.join(settingup.DB_FOLDER, dtype+'_'+settingup.DBXML['Mammalia']), 'wb', 9) as mmoutfile, gzip.open(os.path.join(settingup.DB_FOLDER, dtype+'_'+settingup.DBXML['Fungi']), 'wb', 9) as fgoutfile, gzip.open(os.path.join(settingup.DB_FOLDER, dtype+'_'+settingup.DBXML['Eukaryota']), 'wb', 9) as ekoutfile:
        found_entry = False # Default for entry checkpoint as false.
        
        relative_pos_name = False # Default for checkpoint for entry name.
        sub_loc = False # Default checkpoint for subcellular location.
        recomname_loc = False # Default checkpoint for location.
        goid_loc = False # Default checkpoint for goid.
        feature_loc = False # Default checkpoint for feature.
        feature_loc_seq = False # Default checkpoint for feature location.
        dbr_id = False # Default checkpoint for dbreference.
        ns = '{http://uniprot.org/uniprot}' # Namespace.
        entry = etree.Element('entry') # Create entry xml element.
        protein = ''
        recommendedn = ''
        subcellular_loc = ''
        
        dbreference = ''
        feature_inf = ''
        position_seq = ''
        found_organism = False
        
        outfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        outfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
        
        vboutfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        vboutfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
            
        mmoutfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        mmoutfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
        
        fgoutfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        fgoutfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
        
        ekoutfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        ekoutfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
            
        for event, elem in etree.iterparse(gzip.open(filepath, 'rb'), events = ("start","end")):
            #print(event, etree.tostring(elem))
            #if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and event == 'start':
                #print(etree.tostring(elem))
            
            if elem.tag == ns+'entry' and event == "start":
                found_entry = True 
                
            if elem.tag == ns+'accession' and found_entry == True and event == "end":
                
                accession = etree.SubElement(entry, 'accession')
                accession.text = elem.text
                #print(elem.text)
                #print(etree.tostring(entry))
                relative_pos_name = True
            if elem.tag == ns+'name' and found_entry == True and relative_pos_name == True and event == "end":
                e_name = etree.SubElement(entry,'name')
                e_name.text = elem.text
                
                relative_pos_name = False
            if elem.tag == ns+'protein' and found_entry == True and event == "start":
                protein = etree.SubElement(entry,'protein')
            if elem.tag == ns+'recommendedName' and found_entry == True and event == "start":
                recomname_loc = True
                recommendedn = etree.SubElement(protein, 'recommendedName')
            if elem.tag == ns+'fullName' and recomname_loc == True and event == "end":
                fullname = etree.SubElement(recommendedn, 'fullName')
                fullname.text = elem.text
            if elem.tag == ns+'ecNumber' and recomname_loc == True and event == "end":
                ecnm = etree.SubElement(recommendedn, 'ecNumber')
                ecnm.text = elem.text
            if elem.tag == ns+'recommendedName' and found_entry == True and event == "end":
                recomname_loc = False
            if elem.tag == ns+'organism' and found_entry == True and event == "start":
                found_organism = True
                organism = etree.SubElement(entry, 'organism')
            if elem.tag == ns+'name' and found_entry == True and elem.get('type') == 'scientific' and found_organism == True and event == "end":
                org_name = etree.SubElement(organism, 'name', {'type':'scientific'})
                org_name.text = elem.text
            if elem.tag == ns+'taxon' and found_organism == True and found_entry == True and event == "start" and elem.text in taxa_dict:
                taxa_dict[elem.text] = True
            if elem.tag == ns+'organism' and found_entry == True and event == "end":
                found_organism = False
            
            if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "start":
                sub_loc = True
                subcellular_loc = etree.SubElement(entry, 'comment', {'type': 'subcellular location'})
            
            if elem.tag == ns+'subcellularLocation' and found_entry == True and event == "start":
                
                subr_loc = etree.SubElement(subcellular_loc, 'subcellularLocation')
                  
            if elem.tag == ns+'location' and found_entry == True and sub_loc == True and event == "end":
                location = etree.SubElement(subr_loc, 'location')
                location.text = elem.text
                #print(elem.text)
            if elem.tag == ns+'topology' and found_entry == True and sub_loc == True and event == "end":
                topology = etree.SubElement(subr_loc, 'topology')
                topology.text = elem.text
            if elem.tag == ns+'orientation' and found_entry == True and sub_loc == True and event == "end":
                
                topology = etree.SubElement(subr_loc, 'orientation')
                topology.text = elem.text
                
            if elem.tag == ns+'text' and found_entry == True and sub_loc == True and event == "end":
                text = etree.SubElement(subcellular_loc, 'text')
                text.text = elem.text
                
                #print(elem.text)
            #if elem.tag == ns+'subcellular location' and found_entry == True and event == "end":
                #subcellular_loc.append(subr_loc)
                #subr_loc = etree.Element('subcellular location')
                
            if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "end":
                sub_loc = False
                
            if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "start":
                goid_loc = True
                dbre_at = {'type':elem.get('type'), 'id':elem.get('id')}
                dbReference = etree.SubElement(entry, 'dbReference', dbre_at)
                
            if elem.tag == ns+'property' and found_entry == True and goid_loc == True and event == "end":
                dbre_at = {'type':elem.get('type'), 'value':elem.get('value')}
                dbrprop = etree.SubElement(dbReference, 'property', dbre_at)
            if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "end":
                goid_loc = False
                
            if elem.tag == ns+'feature' and found_entry == True and event == "start":
                feature_loc = True
                fattribute = {'type':elem.get('type')}
                if elem.get('description'):
                    fattribute['description'] = elem.get('description')
                    
                feature_inf = etree.SubElement(entry, 'feature', fattribute)
                
                
            if elem.tag == ns+'location' and found_entry == True and event == "start":
                position_seq = etree.SubElement(feature_inf, 'location')
                feature_loc_seq = True
            if elem.tag == ns+'begin' and found_entry == True and feature_loc_seq == True and event == "end":
                b_start = dict()
                if not 'position' in elem.attrib:
                    b_start['status'] = 'unknown'
                        
                else:
                    b_start['position'] = elem.get('position')
                
                b_pos = etree.SubElement(position_seq, 'begin', b_start)
                
            if elem.tag == ns+'end' and found_entry == True and feature_loc_seq == True and event == "end":
                b_end = dict()
                if not 'position' in elem.attrib:
                    b_end['status'] = 'unknown'
                        
                else:
                    b_end['position'] = elem.get('position')
                e_pos = etree.SubElement(position_seq, 'end', b_end)
                
            if elem.tag == ns+'position' and found_entry == True and feature_loc_seq == True and event == "end":
                b_direct = dict()
                if not 'position' in elem.attrib:
                    b_direct['status'] = 'unknown'
                        
                else:
                    b_direct['position'] = elem.get('position')
                pos = etree.SubElement(position_seq, 'position', b_direct)
                
            if elem.tag == ns+'location' and found_entry == True and event == "end":
                feature_loc_seq = False
                
            if elem.tag == ns+'feature' and found_entry == True and event == "end":
                feature_loc = False
            if elem.tag == ns+'sequence' and found_entry == True and event == "end" and 'checksum' in elem.attrib:
                sattribute = dict()
                if elem.get('length'):
                    sattribute['length'] = elem.get('length')
                if elem.get('mass'):
                    sattribute['mass'] = elem.get('mass')
                sequence = etree.SubElement(entry, 'sequence', sattribute)
                
                sequence.text = elem.text
            if elem.tag == ns+'entry' and found_entry == True and event == "end":
                found_entry = False
                
                #print(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                outfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                if taxa_dict['Vertebrata'] == True:
                    vboutfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                    taxa_dict['Vertebrata'] = False
                if taxa_dict['Mammalia'] == True:
                    mmoutfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                    taxa_dict['Mammalia'] = False
                if taxa_dict['Fungi'] == True:
                    fgoutfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                    taxa_dict['Fungi'] = False
                if taxa_dict['Eukaryota'] == True:
                    ekoutfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                    taxa_dict['Eukaryota'] = False
                count += 1
                print(count, end='\r')
                elem.clear()
                entry = etree.Element('entry')
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
            if elem.tag == ns+'entry' and found_entry == False and event == "end":
                found_entry = False
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
        outfile.write(b'</uniprot>')
        vboutfile.write(b'</uniprot>')
        mmoutfile.write(b'</uniprot>')
        fgoutfile.write(b'</uniprot>')
        ekoutfile.write(b'</uniprot>')
def get_entry(elem, event, ns):
    found_entry = False
    if elem.tag == ns+'entry' and event == "start":
        # Start entry checkpoint.
        found_entry = True
    return found_entry
        
def get_seq(elem, event, ns, **kwargs):
    if elem.tag == ns+'sequence' and event == "end" and 'checksum' in elem.attrib:
        return elem.text
        
def get_acc(elem, event, ns, **kwargs):
    found_entry = kwargs['found_entry']
    if elem.tag == ns+'accession' and event == "end":
        #relative_pos_name = True
        # Start entry name checkpoint.
        # Add each found accession into the temporary entry dictionary.
        return elem.text
    
def get_entryname(elem, event, ns, **kwargs):
    found_entry = kwargs['found_entry']
    relative_pos_name = kwargs['relative_pos_name']
    if elem.tag == ns+'name' and event == "end": 
        #relative_pos_name = False
        return elem.text

def libgen_peptide(seq, rule='trypsin', missc=0):
    results = []
    if seq is not None:
        digested_peps = parser.cleave(seq, parser.expasy_rules[rule], missc)
        for pep in digested_peps:
            result = [pep, mass.calculate_mass(sequence=pep)]
            results.append(result)
        return results
    
#def libgen_pep_db(filepath='uniprot_sprot.xml.gz', dtype='sp'):
    #found_entry = False
    #relative_pos_name = False
    #seq = []
    #pep = []
    #entry_name = []
    #ns = '{http://uniprot.org/uniprot}'
    #with gzip.open(filepath, 'rb') as inputfile:
        #for event, elem in etree.iterparse(inputfile, events = ("start","end")):
            #if found_entry == False:
                #found_entry = get_entry(elem, event, ns)
            #if found_entry == True:
                #acc = get_acc(elem, event, ns)
                #if acc is not None:
            
    #sqlcred = mysql_cred()
    #dbcon = MySQLdb.connect(user=sqlcred['user'], passwd=sqlcred['pass'], host=sqlcred['host'], db=sqlcred['db'])
    #cursor = dbcon.cursor()
    #"""DROP TABLE IF EXISTS """
